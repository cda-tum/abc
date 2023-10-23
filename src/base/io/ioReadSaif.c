/**CFile****************************************************************

  FileName    [ioReadSaif.c]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [Network and node package.]

  Synopsis    [Computes switching activity of nodes in the ABC network.]

  Author      [Benjamin Hien]

  Affiliation [UC Berkeley]

  Date        [Ver. 1.0. Started - June 20, 2005.]

  Revision    [$Id: simSwitch.c,v 1.00 2005/06/20 00:00:00 alanmi Exp $]

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
    int                  t0;           // the time spent at logic '0'
    int                  t1;           // the time spent at logic '1'
    int                  tx;           // the time spent at logic 'X'
    int                  tc;           // the number of total transitions ( '0'->'1' or '1'->'0' )
    int                  ig;           // the number of "instantaneous glitches"
};

typedef struct Saif_SimInfo_t_        Saif_SimInfo_t;   // all reading info
struct Saif_SimInfo_t_
{
    // general simulation info
    int                  Timescale;    // the timescale used (e.g. fs)
    int                  Duration;     // the simulation duration in timescale units
    int                  SimDur;       // the simulation duration in s ( timescale * duration )

    // network instances
    Vec_Ptr_t *          vInstances;   // the line currently parsed
    Vec_Ptr_t *          vInNames;     // the names of the instances
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
    int             nChildren;    // number of children
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

static inline int           Saif_ItemId( Saif_Tree_t * p, Saif_Item_t * pItem )                { return pItem - p->pItems;                                         }

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

static inline int           Saif_Compare( Saif_Tree_t * p, Saif_Pair_t Pair, char * pStr )     { return strncmp( p->pContents+Pair.Beg, pStr, Pair.End-Pair.Beg ) || ((int)strlen(pStr) != Pair.End-Pair.Beg); }

/**Function*************************************************************

  Synopsis    [Returns free item.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
/*int Saif_BuildItem2( Saif_Tree_t * p, char ** ppPos, char * pEnd , int instance_field )
{
    Saif_Item_t * pItem;
    Saif_Pair_t Key, Head, Body;
    char * pNext, * pPrev, * pStop;
    // Skip spaces at the beginning
    ////////////////////////////////////////// Keys: Begin /////////////////////////////////////////
    Key.End = 0;
    if ( Saif_SkipSpaces( p, ppPos, pEnd, 0 ) )
        return -2;
    Key.Beg = *ppPos - p->pContents + 1; // every first character is a '('
    if ( Saif_SkipEntry( ppPos, pEnd ) )
        goto exit;
    Key.End = *ppPos - p->pContents;

    char* test_string = Saif_ReadString( p, Key);
    printf( "Key " );
    printf( test_string );
    printf( "\n" );
    ////////////////////////////////////////// Keys: End ///////////////////////////////////////////

    pNext = *ppPos;
    pPrev = pNext - 1;

    if (*pNext == '\n')
    {

        pItem = Saif_NewItem(p);
        pItem->Key = Key;
        //pItem->Head = Saif_UpdateHead( p, Head );
        printf("Child ");
        printf("\n");
        pItem->Child = Saif_BuildItem2(p, ppPos, pEnd, 0);
        if (pItem->Child == -1)
            goto exit;
        return Saif_ItemId(p, pItem);
    }

    if (*pNext == ' ')
    {
        Saif_SkipSpaces( p, ppPos, pEnd, 0 );
        pNext = *ppPos;
        if (*pNext == ')')
        {
            (*ppPos) += 1;
            pItem = Saif_NewItem(p);
            pItem->Key = Key;
            printf("Next1 ");
            printf("\n");
            //pItem->Head = Saif_UpdateHead( p, Head );
            pItem->Next = Saif_BuildItem2(p, ppPos, pEnd, 0);
            if (pItem->Next == -1)
                goto exit;
            return Saif_ItemId(p, pItem);
        }
        if (*pNext == '\"')
        {
            Head.Beg = *ppPos - p->pContents;
            if (Saif_SkipEntry(ppPos, pEnd))
                goto exit;
            Head.End = *ppPos - p->pContents;
            *//*if ( Saif_SkipSpaces( p, ppPos, pEnd, 0 ) )
                goto exit;*//*
            (*ppPos) += 1;

            char * test_string = Saif_ReadString(p, Head);

            printf("Head ");
            printf(test_string);
            printf("\n");
            pItem = Saif_NewItem(p);
            printf("Next2 ");
            printf("\n");
            pItem->Key = Key;
            pItem->Head = Saif_UpdateHead(p, Head);
            pItem->Next = Saif_BuildItem2(p, ppPos, pEnd, 0 );
            return Saif_ItemId(p, pItem);
        }

        Head.Beg = *ppPos - p->pContents;
        if (Saif_SkipEntry(ppPos, pEnd))
            goto exit;
        Head.End = *ppPos - p->pContents;
        pNext = *ppPos;
        pPrev = pNext - 1;
        if (*pNext == ' ')
        {
            Saif_SkipSpaces(p, ppPos, pEnd, 0);
            pNext = *ppPos;
            if (*pNext == ')')
            {
                char * test_string = Saif_ReadString(p, Head);
                printf("Head ");
                printf(test_string);
                printf("\n");
                (*ppPos) += 1;
                pItem = Saif_NewItem(p);
                printf("Next3 ");
                printf("\n");
                pItem->Key = Key;
                pItem->Head = Saif_UpdateHead( p, Head );
                pItem->Next = Saif_BuildItem2(p, ppPos, pEnd, 0);
                if (pItem->Next == -1)
                    goto exit;
                return Saif_ItemId(p, pItem);
            }
        }
        // Also check previous for being a ')': If so were in the T field
        if (*pNext == '\n' || *pNext == '(')
        {
            if (*pNext == '\n' && *pPrev != ')')
            {
                printf("HERE! ");
                printf("\n");
            }
            Head.End -= 1;
            char * test_string = Saif_ReadString(p, Head);
            printf("Head ");
            printf(test_string);
            printf("\n");
            pItem = Saif_NewItem(p);
            printf("Next4 ");
            printf("\n");
            pItem->Key = Key;
            pItem->Head = Saif_UpdateHead(p, Head);
            pItem->Next = Saif_BuildItem2(p, ppPos, pEnd, 0);
            if (pItem->Next == -1)
                goto exit;
            return Saif_ItemId(p, pItem);
        }
        Body.Beg = *ppPos - p->pContents;
        Saif_SkipEntry( ppPos, pEnd );
        Body.End = *ppPos - p->pContents - 1;
        pNext = *ppPos;
        char * test_string = Saif_ReadString(p, Head);
        printf("Head ");
        printf(test_string);
        printf("\n");
        test_string = Saif_ReadString(p, Body);
        printf("Body ");
        printf(test_string);
        printf("\n");
        pItem = Saif_NewItem(p);
        pItem->Key = Key;
        printf("Next5 ");
        printf("\n");
        pItem->Head = Saif_UpdateHead(p, Head);
        pItem->Body = Body;
        pItem->Next = Saif_BuildItem2(p, ppPos, pEnd, 0);
        if (pItem->Next == -1)
            goto exit;
        return Saif_ItemId(p, pItem);
    }

    exit:
    if ( p->pError == NULL )
    {
        p->pError = ABC_ALLOC( char, 1000 );
        sprintf( p->pError, "File \"%s\". Line %6d. Failed to parse entry \"%s\".\n",
                 p->pFileName, p->nLines, Saif_ReadString(p, Key) );
    }
    return -1;
}*/

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
void Saif_StringDump( char * pFileName, Vec_Str_t * vStr )
{
    FILE * pFile = fopen( pFileName, "wb" );
    int RetValue = 0;
    if ( pFile == NULL )
    {
        printf( "Saif_StringDump(): The output file is unavailable.\n" );
        return;
    }
    RetValue = fwrite( Vec_StrArray(vStr), 1, Vec_StrSize(vStr), pFile );
    fclose( pFile );
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

Vec_Int_t * Io_ReadSaif( char * pFileName )
{
    Vec_Int_t * vSwitching;
    float * pSwitching;
    Vec_Ptr_t * vNodes;
    Vec_Ptr_t * vSimInfo;
    Abc_Obj_t * pNode;
    unsigned * pSimInfo;
    int nSimWords, i;

    Saif_Tree_t * p;
    char * pPos;

    // A network has to be read im
    //pPos = p->pContents;
    p = Saif_Parse( pFileName );
    if ( p == NULL )
        return NULL;
    Saif_Stop( p, 0 );

    // This function creates the
    // Saif_CreateInfo

    printf("Function executed");

    // The network has to correspond to the saif data

    return 0;
}