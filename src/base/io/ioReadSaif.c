/**CFile****************************************************************

  FileName    [ioReadSaif.c]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [Network and node package.]

  Synopsis    [Returns switching activity of nodes in the ABC network from SAIF file.]

  Author      [Benjamin Hien]

  Affiliation [TU Munich]

  Date        [Ver. 1.0. Started - Oct 17, 2023.]

  Revision    []

***********************************************************************/

#include "base/abc/abc.h"

ABC_NAMESPACE_IMPL_START


////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

typedef struct Saif_Pin_t_        Saif_Pin_t;   // all reading info
struct Saif_Pin_t_
{
    // switching information
    char *               pin_name;     // pin name
    long                 t0;           // the time spent at logic '0'
    long                 t1;           // the time spent at logic '1'
    long                 tx;           // the time spent at logic 'X'
    long                 tc;           // the number of total transitions ( '0'->'1' or '1'->'0' )
    long                 ig;           // the number of "instantaneous glitches"
};

typedef struct Saif_Instance_t_        Saif_Instance_t;   // all reading info
struct Saif_Instance_t_
{
    // general simulation info
    char *               Name;         // the timescale used (e.g. fs)

    // network instances
    Vec_Ptr_t            vPins;           // NamedSet<Saif_Pin_t>
    int                  n_inputs;        // -- 'pins[0 … n_inputs-1]' are input pins
    int                  n_outputs;       // -- 'pins[n_inputs … n_inputs+n_outputs-1]' are output pins
    int                  Id;              // instance ID
    int                  Top_Id;          // connected instance  higher in hierarchy
    Vec_Int_t            Sub_Id;          // connected instances lower in hierarchy
};

typedef struct Saif_SimInfo_t_        Saif_SimInfo_t;   // all reading info
struct Saif_SimInfo_t_
{
    // general simulation info
    int                  Timescale;    // the timescale used (e.g. fs)
    long                 Duration;     // the simulation duration in timescale units
    int                  SimDur;       // the simulation duration in s ( timescale * duration )

    // network instances
    Vec_Ptr_t            vInstances;   // NamedSet<Saif_Instance_t>
};

typedef struct Saif_Pair_t_ Saif_Pair_t;
struct Saif_Pair_t_
{
    int             Beg;          // item beginning
    int             End;          // item end
};

typedef struct Saif_Item_t_ Saif_Item_t;
struct Saif_Item_t_
{
    int             iLine;        // file line where the item's spec begins
    Saif_Pair_t     Key;          // key part
    Saif_Pair_t     Head;         // head part
    Saif_Pair_t     Body;         // body part
    int             Next;         // next item in the list
    int             Child;        // first child item
};

typedef struct Saif_Tree_t_ Saif_Tree_t;
struct Saif_Tree_t_
{
    char *          pFileName;    // input Liberty file name
    char *          pContents;    // file contents
    int             nContents;    // file size
    int             nLines;       // line counter
    int             nItems;       // number of items
    int             nItermAlloc;  // number of items allocated
    Saif_Item_t *   pItems;       // the items
    char *          pError;       // the error string
    abctime         clkStart;     // beginning time
    Vec_Str_t *     vBuffer;      // temp string buffer
};

typedef struct Io_ReadSaif_t_        Io_ReadSaif_t;   // all reading info
struct Io_ReadSaif_t_
{
    // general info about file
    char *               pFileName;    // the name of the file
    char *               pVersion;     // the version of SAIF
    // current processing info
    int                  LineCur;      // the line currently parsed
    // temporary storage for tokens
    Vec_Ptr_t *          vTokens;      // the current tokens
    Vec_Ptr_t *          vNewTokens;   // the temporary storage for the tokens
    Vec_Str_t *          vCubes;       // the temporary storage for the tokens

    // the error message
    FILE *               Output;       // the output stream
    char                 sError[1000]; // the error string generated during parsing
    int                  fError;       // set to 1 when error occurs
};
////////////////////////////////////////////////////////////////////////
///                     FUNCTION DEFINITIONS                         ///
////////////////////////////////////////////////////////////////////////


/**Function*************************************************************

  Synopsis    [Gets the name to write.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

static inline Saif_Item_t *  Saif_Root( Saif_Tree_t * p )                                      { return p->pItems;                                                 }
static inline Saif_Item_t *  Saif_Item( Saif_Tree_t * p, int v )                               { assert( v < p->nItems ); return v < 0 ? NULL : p->pItems + v;     }
static inline int            Saif_Compare( Saif_Tree_t * p, Saif_Pair_t Pair, char * pStr )    { return strncmp( p->pContents+Pair.Beg, pStr, Pair.End-Pair.Beg ) || ((int)strlen(pStr) != Pair.End-Pair.Beg); }
static inline void           Saif_PrintWord( FILE * pFile, Saif_Tree_t * p, Saif_Pair_t Pair )        { char * pBeg = p->pContents+Pair.Beg, * pEnd = p->pContents+Pair.End; while ( pBeg < pEnd ) fputc( *pBeg++, pFile ); }
static inline void           Saif_PrintSpace( FILE * pFile, int nOffset )                             { int i; for ( i = 0; i < nOffset; i++ ) fputc(' ', pFile);         }
static inline int            Saif_ItemId( Saif_Tree_t * p, Saif_Item_t * pItem )                     { return pItem - p->pItems;                                         }

#define Saif_ItemForEachChild( p, pItem, pChild ) \
    for ( pChild = Saif_Item(p, pItem->Child); pChild; pChild = Saif_Item(p, pChild->Next) )
#define Saif_ItemForEachChildName( p, pItem, pChild, pName ) \
    for ( pChild = Saif_Item(p, pItem->Child); pChild; pChild = Saif_Item(p, pChild->Next) ) if ( Saif_Compare(p, pChild->Key, pName) ) {} else

int Saif_CountItems( char * pBeg, char * pEnd )
{
    int Counter = 0;
    for ( ; pBeg < pEnd; pBeg++ )
        Counter += (*pBeg == '(' );
    return Counter;
}
// removes C-style comments
/*
void Saif_WipeOutComments( char * pBeg, char * pEnd )
{
    char * pCur, * pStart;
    for ( pCur = pBeg; pCur < pEnd; pCur++ )
    if ( pCur[0] == '/' && pCur[1] == '*' )
        for ( pStart = pCur; pCur < pEnd; pCur++ )
        if ( pCur[0] == '*' && pCur[1] == '/' )
        {
            for ( ; pStart < pCur + 2; pStart++ )
            if ( *pStart != '\n' ) *pStart = ' ';
            break;
        }
}
*/
void Saif_WipeOutComments( char * pBeg, char * pEnd )
{
    char * pCur, * pStart;
    for ( pCur = pBeg; pCur < pEnd-1; pCur++ )
        if ( pCur[0] == '/' && pCur[1] == '*' )
        {
            for ( pStart = pCur; pCur < pEnd-1; pCur++ )
                if ( pCur[0] == '*' && pCur[1] == '/' )
                {
                    for ( ; pStart < pCur + 2; pStart++ )
                        if ( *pStart != '\n' ) *pStart = ' ';
                    break;
                }
        }
        else if ( pCur[0] == '/' && pCur[1] == '/' )
        {
            for ( pStart = pCur; pCur < pEnd; pCur++ )
                if ( pCur[0] == '\n' || pCur == pEnd-1 )
                {
                    for ( ; pStart < pCur; pStart++ ) *pStart = ' ';
                    break;
                }
        }
}
static inline int Saif_CharIsSpace( char c )
{
    return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\\';
}
static inline int Saif_SkipSpaces( Saif_Tree_t * p, char ** ppPos, char * pEnd, int fStopAtNewLine )
{
    char * pPos = *ppPos;
    for ( ; pPos < pEnd; pPos++ )
    {
        if ( *pPos == '\n' )
        {
            p->nLines++;
            if ( fStopAtNewLine )
                break;
        }
        if ( !Saif_CharIsSpace(*pPos) )
            break;
    }
    *ppPos = pPos;
    return pPos == pEnd;
}
// skips entry delimited by " :;(){}" and returns 1 if reached the end
static inline int Saif_SkipEntry( char ** ppPos, char * pEnd )
{
    char * pPos = *ppPos;
    if ( *pPos == '\"' )
    {
        for ( pPos++; pPos < pEnd; pPos++ )
            if ( *pPos == '\"' )
            {
                pPos++;
                break;
            }
    }
    else
    {
        for ( ; pPos < pEnd; pPos++ )
            if ( *pPos == ' ' || *pPos == '\r' || *pPos == '\n' || *pPos == '\t' ||
                 *pPos == '(' || *pPos == ')' )
                break;
    }
    *ppPos = pPos;
    return pPos == pEnd;
}
// finds the matching closing symbol
static inline char * Saif_FindMatch( char * pPos, char * pEnd )
{
    int Counter = 1;

    if ( *pPos == '\n' || *pPos == '\"' || *pPos == ' ' )
    {
        for ( ; pPos < pEnd; pPos++ )
        {
            if ( *pPos == '(' )
                Counter++;
            if ( *pPos == ')' )
                Counter--;
            if ( Counter == 0 )
                break;
        }
    }
    else if ( *pPos == ')' )
    {
        Counter = 0;
        for ( ; pPos < pEnd; pPos++ )
        {
            if ( *pPos == '(' )
                Counter++;
            if ( *pPos == ')' )
                Counter--;
            if ( Counter == 0 )
                break;
        }
    }


    assert( *pPos == ')' );
    return pPos;
}
// trims spaces around the head
static inline Saif_Pair_t Saif_UpdateHead( Saif_Tree_t * p, Saif_Pair_t Head )
{
    Saif_Pair_t Res;
    char * pBeg = p->pContents + Head.Beg;
    char * pEnd = p->pContents + Head.End;
    char * pFirstNonSpace = NULL;
    char * pLastNonSpace = NULL;
    char * pChar;
    for ( pChar = pBeg; pChar < pEnd; pChar++ )
    {
        if ( *pChar == '\n' )
            p->nLines++;
        if ( Saif_CharIsSpace(*pChar) )
            continue;
        pLastNonSpace = pChar;
        if ( pFirstNonSpace == NULL )
            pFirstNonSpace = pChar;
    }
    if ( pFirstNonSpace == NULL || pLastNonSpace == NULL )
        return Head;
    assert( pFirstNonSpace && pLastNonSpace );
    Res.Beg = pFirstNonSpace - p->pContents;
    Res.End = pLastNonSpace  - p->pContents + 1;
    return Res;
}

static inline Saif_Item_t * Saif_NewItem( Saif_Tree_t * p )
{
    p->pItems[p->nItems].iLine = p->nLines;
    p->pItems[p->nItems].Child = -1;
    p->pItems[p->nItems].Next  = -1;
    return p->pItems + p->nItems++;
}

/**Function*************************************************************

  Synopsis    [Gets the name to write.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
char * Saif_ReadString( Saif_Tree_t * p, Saif_Pair_t Pair )
{
    // static char Buffer[ABC_MAX_LIB_STR_LEN];
    char * Buffer;
    if ( Pair.End - Pair.Beg + 2 > Vec_StrSize(p->vBuffer) )
        Vec_StrFill( p->vBuffer, Pair.End - Pair.Beg + 100, '\0' );
    Buffer = Vec_StrArray( p->vBuffer );
    strncpy( Buffer, p->pContents+Pair.Beg, Pair.End-Pair.Beg );
    if ( Pair.Beg < Pair.End && Buffer[0] == '\"' )
    {
        assert( Buffer[Pair.End-Pair.Beg-1] == '\"' );
        Buffer[Pair.End-Pair.Beg-1] = 0;
        return Buffer + 1;
    }
    Buffer[Pair.End-Pair.Beg] = 0;
    return Buffer;
}
/**Function*************************************************************

  Synopsis    [Returns free item.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int Saif_BuildItem2(Saif_Tree_t * p, char ** ppPos, char * pEnd ) {
    Saif_Item_t * pItem;
    Saif_Pair_t Key, Body;
    char * pNext, * pStop;

    if (Saif_SkipSpaces(p, ppPos, pEnd, 0) == 1)
        return -2;

    Key.Beg = *ppPos - p->pContents;
    if (Saif_SkipEntry(ppPos, pEnd) == 1)
        goto exit;

    Key.End = *ppPos - p->pContents;
    pStop = Saif_FindMatch(*ppPos, pEnd);

    Body.Beg = *ppPos - p->pContents + 1;
    Body.End = pStop - p->pContents;

    switch (**ppPos) {
        case '\n':
            pItem = Saif_NewItem(p);
            pItem->Key = Key;
            pItem->Body = Body;
            *ppPos += 1;
            pItem->Child = Saif_BuildItem2(p, ppPos, pStop);
            if (pItem->Child == -1) goto exit;
            *ppPos = pStop + 1;
            pItem->Next = Saif_BuildItem2(p, ppPos, pEnd);
            if (pItem->Next == -1) goto exit;
            return Saif_ItemId(p, pItem);

        case ')':
            *ppPos += 2;
            pItem = Saif_NewItem(p);
            pItem->Key = Key;
            pItem->Next = Saif_BuildItem2(p, ppPos, pEnd);
            if (pItem->Next == -1) goto exit;
            return Saif_ItemId(p, pItem);

        case '\"':
            if (Saif_SkipEntry(ppPos, pEnd)) goto exit;
            // Continue to default.

        default:
            *ppPos += 1;  // Skip the current character.
            if (**ppPos == '\n') {
                pItem = Saif_NewItem(p);
                Saif_Pair_t Head;
                Head.Beg = *ppPos - p->pContents;
                if (Saif_SkipEntry(ppPos, pEnd)) goto exit;
                Head.End = *ppPos - p->pContents;
                pItem->Head = Saif_UpdateHead(p, Head);
                pItem->Key = Key;
                pItem->Body = Body;
                pItem->Child = Saif_BuildItem2(p, ppPos, pStop);
                if (pItem->Child == -1) goto exit;
                *ppPos = pStop + 1;
                pItem->Next = Saif_BuildItem2(p, ppPos, pEnd);
                if (pItem->Next == -1) goto exit;
                return Saif_ItemId(p, pItem);
            }
            if (Saif_SkipSpaces(p, ppPos, pEnd, 0))
                goto exit;

            pItem = Saif_NewItem(p);
            pItem->Key  = Key;
            pItem->Body = Body;
            pItem->Next = Saif_BuildItem2( p, ppPos, pEnd );
            if (pItem->Next == -1)
                goto exit;
            return Saif_ItemId(p, pItem);
    }

    exit:
    if ( p->pError == NULL )
    {
        p->pError = ABC_ALLOC( char, 1000 );
        sprintf( p->pError, "File \"%s\". Line %6d. Failed to parse entry \"%s\".\n",
                 p->pFileName, p->nLines, Saif_ReadString( p, Key) );
    }
    return -1;
}

/**Function*************************************************************

  Synopsis    [Returns free item.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int Saif_BuildItem( Saif_Tree_t * p, char ** ppPos, char * pEnd )
{
    Saif_Item_t * pItem;
    Saif_Pair_t Key, Head, Body;
    char * pNext, * pStop;
    if ( Saif_SkipSpaces( p, ppPos, pEnd, 0 ) )
        return -2;
    pNext = *ppPos;
    assert( *pNext == '(' );
    if ( *pNext == '(' )
    {
        (*ppPos) += 1;
    }
    Key.End = 0;
    Key.Beg = *ppPos - p->pContents; // +1 first char is always '('
    if ( Saif_SkipEntry( ppPos, pEnd ) )
        goto exit;
    Key.End = *ppPos - p->pContents;
    /*if ( Saif_SkipSpaces( p, ppPos, pEnd, 0 ) )
        goto exit;*/
    pNext = *ppPos;
    // end of body
    // pNext = *ppPos;
    pStop = Saif_FindMatch( pNext, pEnd );
    Body.Beg = pNext - p->pContents + 1;
    Body.End = pStop - p->pContents;
    // First Entry, NET, Pins
    if ( *pNext == '\n' )
    {
        pItem = Saif_NewItem( p );
        pItem->Key  = Key;
        // pItem->Head = Saif_UpdateHead( p, Head );
        pItem->Body = Body;
        *ppPos = pNext + 1;
        pItem->Child = Saif_BuildItem( p, ppPos, pStop );
        if ( pItem->Child == -1 )
            goto exit;
        *ppPos = pStop + 1;
        pItem->Next = Saif_BuildItem( p, ppPos, pEnd );
        if ( pItem->Next == -1 )
            goto exit;
        return Saif_ItemId( p, pItem );
    }
    if ( Saif_SkipSpaces( p, ppPos, pEnd, 0 ) )
        goto exit;
    pNext = *ppPos;
    // Only Design
    if ( *pNext == ')' )
    {
        *ppPos = pNext + 2;
        pItem = Saif_NewItem( p );
        pItem->Key  = Key;
        pItem->Next = Saif_BuildItem( p, ppPos, pEnd );
        if ( pItem->Next == -1 )
            goto exit;
        return Saif_ItemId( p, pItem );
    }
    if (*pNext == '\"')
    {
        if ( Saif_SkipEntry( ppPos, pEnd ) )
            goto exit;
    }
    Head.Beg = *ppPos - p->pContents;
    if ( Saif_SkipEntry( ppPos, pEnd ) )
        goto exit;
    Head.End = *ppPos - p->pContents;
    pNext = *ppPos;
    // For INSTANCE
    if ( *pNext == '\n' )
    {
        pItem = Saif_NewItem( p );
        pItem->Key  = Key;
        pItem->Head = Saif_UpdateHead( p, Head );
        pItem->Body = Body;
        *ppPos = pNext + 1;
        pItem->Child = Saif_BuildItem( p, ppPos, pStop );
        if ( pItem->Child == -1 )
            goto exit;
        *ppPos = pStop + 1;
        pItem->Next = Saif_BuildItem( p, ppPos, pEnd );
        if ( pItem->Next == -1 )
            goto exit;
        return Saif_ItemId( p, pItem );
    }
    // end of list
    if ( Saif_SkipSpaces( p, ppPos, pEnd, 1 ) )
        goto exit;
    pNext = *ppPos;
    if ( *pNext == ')' )
        *ppPos = pNext + 2;
    else
    {
        if ( Saif_SkipEntry( ppPos, pEnd ) )
            goto exit;
        // end of list
        if ( Saif_SkipSpaces( p, ppPos, pEnd, 1 ) )
            goto exit;
        pNext = *ppPos;
    }
    if ( *pNext == ')' )
        *ppPos = pNext + 2;
    // For last entry of T
    if ( Saif_SkipSpaces( p, ppPos, pEnd, 0 ) )
    {
        pItem = Saif_NewItem( p );
        pItem->Key  = Key;
        pItem->Body = Body;
        return Saif_ItemId( p, pItem );
    }
    if ( Saif_SkipSpaces( p, ppPos, pEnd, 0 ) )
        goto exit;
    pItem = Saif_NewItem( p );
    pItem->Key  = Key;
    pItem->Body = Body;
    pItem->Next = Saif_BuildItem( p, ppPos, pEnd );
    if ( pItem->Next == -1 )
        goto exit;
    return Saif_ItemId( p, pItem );

    exit:
    if ( p->pError == NULL )
    {
        p->pError = ABC_ALLOC( char, 1000 );
        sprintf( p->pError, "File \"%s\". Line %6d. Failed to parse entry \"%s\".\n",
                 p->pFileName, p->nLines, Saif_ReadString(p, Key) );
    }
    return -1;
}

/**Function*************************************************************

  Synopsis    [File management.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Saif_FixFileName( char * pFileName )
{
    char * pHead;
    for ( pHead = pFileName; *pHead; pHead++ )
        if ( *pHead == '>' )
            *pHead = '\\';
}
int Saif_FileSize( char * pFileName )
{
    FILE * pFile;
    int nFileSize;
    pFile = fopen( pFileName, "rb" );
    if ( pFile == NULL )
    {
        printf( "Saif_FileSize(): The input file is unavailable (absent or open).\n" );
        return 0;
    }
    fseek( pFile, 0, SEEK_END );
    nFileSize = ftell( pFile );
    fclose( pFile );
    return nFileSize;
}
char * Saif_FileContents( char * pFileName, int nContents )
{
    FILE * pFile = fopen( pFileName, "rb" );
    char * pContents = ABC_ALLOC( char, nContents+1 );
    int RetValue = 0;
    RetValue = fread( pContents, nContents, 1, pFile );
    fclose( pFile );
    pContents[nContents] = 0;
    return pContents;
}

/**Function*************************************************************

  Synopsis    [Starts the parsing manager.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
Saif_Tree_t * Saif_Start( char * pFileName )
{
    Saif_Tree_t * p;
    int RetValue;
    // read the file into the buffer
    Saif_FixFileName( pFileName );
    RetValue = Saif_FileSize( pFileName );
    if ( RetValue == 0 )
        return NULL;
    // start the manager
    p = ABC_ALLOC( Saif_Tree_t, 1 );
    memset( p, 0, sizeof(Saif_Tree_t) );
    p->clkStart  = Abc_Clock();
    p->nContents = RetValue;
    p->pContents = Saif_FileContents( pFileName, p->nContents );
    // other
    p->pFileName = Abc_UtilStrsav( pFileName );
    // CountItems adjusted
    p->nItermAlloc = 10 + Saif_CountItems( p->pContents, p->pContents+p->nContents );
    p->pItems = ABC_CALLOC( Saif_Item_t, p->nItermAlloc );
    p->nItems = 0;
    p->nLines = 1;
    p->vBuffer = Vec_StrStart( 10 );
    return p;
}
void Saif_Stop( Saif_Tree_t * p, int fVerbose )
{
    if ( fVerbose )
    {
        printf( "Memory = %7.2f MB. ", 1.0 * (p->nContents + p->nItermAlloc * sizeof(Saif_Item_t))/(1<<20) );
        ABC_PRT( "Time", Abc_Clock() - p->clkStart );
    }
    Vec_StrFree( p->vBuffer );
    ABC_FREE( p->pFileName );
    ABC_FREE( p->pContents );
    ABC_FREE( p->pItems );
    ABC_FREE( p->pError );
    ABC_FREE( p );
}
Saif_Tree_t * Saif_Parse( char * pFileName )
{
    Saif_Tree_t * p;
    char * pPos;
    // Allocate Memory Saif_Tree_t
    if ( (p = Saif_Start(pFileName)) == NULL )
        return NULL;
    pPos = p->pContents;
    // Todo(hien): Check if comments are used in saif files
    Saif_WipeOutComments( p->pContents, p->pContents+p->nContents );
    if ( (!Saif_BuildItem( p, &pPos, p->pContents + p->nContents )) == 0 )
    {
        if ( p->pError ) printf( "%s", p->pError );
        printf( "Parsing failed.  " );
        Abc_PrintTime( 1, "Parsing time", Abc_Clock() - p->clkStart );
    }

    return p;
}

static inline Saif_SimInfo_t * Abc_SaifAlloc()
{
    Saif_SimInfo_t * p;
    p = ABC_CALLOC( Saif_SimInfo_t, 1 );
    p->Timescale = -1;
    p->Duration = -1;
    p->SimDur = -1;
    return p;
}
// This is needed for each vector to dynamically allocate memory
static inline Saif_Instance_t * Abc_SaifInstanceAlloc()
{
    Saif_Instance_t * p;
    p = ABC_CALLOC( Saif_Instance_t, 1 );
    // p->Name = "";
    p->n_inputs = -1;
    p->n_outputs = -1;
    p->Id = -1;
    p->Top_Id = -1;
    return p;
}

static inline Saif_Pin_t * Abc_SaifPinAlloc()
{
    Saif_Pin_t * p;
    p = ABC_CALLOC( Saif_Pin_t, 1 );
    // p->pin_name = "";
    p->t0 = -1;
    p->t1 = -1;
    p->tx = -1;
    p->tc = -1;
    p->ig = -1;
    return p;
}

#define Saif_ForEachInstance( p, pWL, i )       Vec_PtrForEachEntry( Saif_Instance_t *, &p->vInstances, pWL, i )
#define Saif_ForEachPin( p, pWL, i )            Vec_PtrForEachEntry( Saif_Pin_t *, &p->vPins, pWL, i )

static inline void Abc_SaifPinFree( Saif_Pin_t * p )
{
    free(p->pin_name);
    ABC_FREE( p );
}
static inline void Abc_SaifInstanceFree( Saif_Instance_t * p )
{
    int i;
    Saif_Pin_t * pPin;
    free(p->Name);
    // ID Vector needs to be erased
    // p->Sub_Id;
    Vec_IntErase(&p->Sub_Id);
    Saif_ForEachPin( p, pPin, i )
        Abc_SaifPinFree( pPin );
    Vec_PtrErase( &p->vPins );
    ABC_FREE( p );
}
static inline void Abc_SaifFree( Saif_SimInfo_t * p )
{
    int i;
    Saif_Instance_t * pIns;
    Saif_ForEachInstance( p, pIns, i )
        Abc_SaifInstanceFree( pIns );
    Vec_PtrErase( &p->vInstances );
    ABC_FREE( p );
}

// Read out the switching activity information of each pin
void Saif_ReadNetInfo( Saif_Tree_t * p, Saif_Item_t * pRootItem, Saif_Instance_t * pInstance )
{
    Saif_Item_t * pItemNet, * pItemPin, * pItemT;
    // NET items
    Saif_ItemForEachChildName( p, pRootItem, pItemNet, "NET" )
    {
        // Pin items
        Saif_ItemForEachChild( p, pItemNet, pItemPin )
        {
            // Allocate Memory
            Saif_Pin_t * pPin = Abc_SaifPinAlloc();
            Vec_PtrPush( &pInstance->vPins, pPin );
            // Write Name
            pPin->pin_name = Abc_UtilStrsav( Saif_ReadString(p, pItemPin->Key) );
            // Write T(ime) entries
            long* pPinMembers[] = {&pPin->t0, &pPin->t1, &pPin->tx, &pPin->tc, &pPin->ig};
            char* keys[] = {"T0", "T1", "TX", "TC", "IG"};
            int size = sizeof(keys)/sizeof(keys[0]);
            // T(ime) items
            Saif_ItemForEachChild( p, pItemPin, pItemT )
            {
                for(int i=0; i<size; i++)
                {
                    if( !Saif_Compare(p, pItemT->Key, keys[i]) )
                    {
                        *pPinMembers[i] = atof(Saif_ReadString(p, pItemT->Body));
                    }
                }
            }
        }
    }
}
void InstanceHandling(Saif_Instance_t * pInstance, Saif_Item_t * pItem, Vec_Int_t * vTemp, Saif_Tree_t * p, int * iT, int * iC, Saif_SimInfo_t * vSimInf)
{
    // First Instance
    if ( Vec_IntSize(vTemp) == 0 )
    {
        Vec_IntPush( vTemp, *iT );
    }
    else
    {
        pInstance->Top_Id =  Vec_IntEntry( vTemp, *iT );
        if (pItem->Next < 0)
        {
            if (pItem->Child < 0)
            {
                (*iT)--;
            }
            else
            {
                if (Saif_Item(p, pItem->Child)->Next >= 0)
                {
                    (*iT)++;
                    Vec_IntSetEntry( vTemp, *iT, *iC );
                }
                else
                {
                    (*iT)--;
                }
            }
        }
        else
        {
            if (pItem->Child >= 0)
            {
                if (Saif_Item(p, pItem->Child)->Next >= 0)
                {
                    (*iT)++;
                    Vec_IntSetEntry( vTemp, *iT, *iC );
                }
            }
        }
    }
    if( pInstance->Top_Id >= 0 )
    {
        Saif_Instance_t * TopInstance = Vec_PtrEntry( &vSimInf->vInstances, pInstance->Top_Id );
        Vec_IntPush( &TopInstance->Sub_Id, *iC );
    }
    (*iC)++;
}
// Gets as Input the top instance(s)
// Iterates through all  the sub instances and reads out the switching activity information of each
void Saif_ReadInstances( Saif_Tree_t * p, Saif_Item_t * pRootItem, Saif_SimInfo_t * vSimInf, Vec_Int_t * vTemp, int * iT, int * iC )
{
    Saif_Item_t * pItem;

    // Iterate through the Instances
    Saif_ItemForEachChildName( p, pRootItem, pItem, "INSTANCE" )
    {
        Saif_Instance_t * pInstance = Abc_SaifInstanceAlloc();
        Vec_PtrPush( &vSimInf->vInstances, pInstance );
        // Assign the name
        pInstance->Name = Abc_UtilStrsav(Saif_ReadString(p, pItem->Head));
        // Assign the incoming and outgoing instances
        InstanceHandling(pInstance, pItem, vTemp, p, iT, iC, vSimInf);
        // Read in the Net Info (Time information)
        Saif_ReadNetInfo ( p, pItem, pInstance );
        // Recursive call
        Saif_ReadInstances( p, pItem, vSimInf, vTemp, iT, iC );
    }
}

int Saif_ReadTimescale( Saif_Tree_t * p )
{
    Saif_Item_t * pItem;
    Saif_ItemForEachChildName( p, Saif_Root(p), pItem, "TIMESCALE" )
    {
        int unit_zeros[6] = {15, 12,9, 6, 3, 0};
        char *units[6] = {"fs","ps","ns", "us", "ms", "s"};
        char *token;
        long number;
        int number_zeros = 0;
        char *unit = NULL;
        char * timescale = Saif_ReadString(p, pItem->Body);

        // Get the first token (number)
        token = strtok(timescale, " ");
        if (token != NULL) {
            number = strtol(token, NULL, 10);
            while (number != 0) {
                number /= 10;
                number_zeros++;
            }
            number_zeros--;
        }

        // Get the second token (unit)
        token = strtok(NULL, " ");
        if (token != NULL) {
            for(int i = 0; i < 6; i++) {
                if(strcmp(token, units[i]) == 0)
                    return unit_zeros[i] - number_zeros;
            }
        }

        return -1; // Return -1 in case of an unrecognized unit
    }
    printf( "Libery parser cannot read \"time_unit\".  Assuming   time_unit : \"1ns\".\n" );
    return 9;
}

long Saif_ReadDuration( Saif_Tree_t * p )
{
    Saif_Item_t * pItem;
    Saif_ItemForEachChildName( p, Saif_Root(p), pItem, "DURATION" )
        return atol(Saif_ReadString(p, pItem->Body));
    return 10000000000;
}


Saif_SimInfo_t * Saif_CreateInfo( Saif_Tree_t * p )
{
    // Allocate Memory
    Saif_SimInfo_t * vSimInf;
    Vec_Int_t * vTemp;
    int iT = 0, iC = 0;
    vSimInf = Abc_SaifAlloc();
    vTemp = Vec_IntAlloc(10);

    vSimInf->Timescale = Saif_ReadTimescale( p );
    vSimInf->Duration = Saif_ReadDuration( p );
    vSimInf->SimDur = 0; // adjust

    Saif_ReadInstances( p,  Saif_Root(p), vSimInf, vTemp, &iT, &iC );

    Vec_IntFree( vTemp );

    return vSimInf;
}

Vec_Int_t * Io_ReadSaif( char * pFileName )
{
    Saif_Tree_t * p;
    char * pPos;

    Saif_SimInfo_t * vSwitching;

    // A network has to be read im
    //pPos = p->pContents;
    p = Saif_Parse( pFileName );
    if ( p == NULL )
        return NULL;
    vSwitching = Saif_CreateInfo( p );
    // Free Memory Saif_SimInfo_t
    Abc_SaifFree ( vSwitching );
    // Free Memory Saif_Tree_t
    Saif_Stop( p, 0 );

    // This function creates the
    // Saif_CreateInfo

    printf("Function executed");

    // The network has to correspond to the saif data

    return 0;
}

////////////////////////////////////////////////////////////////////////
///                       END OF FILE                                ///
////////////////////////////////////////////////////////////////////////


