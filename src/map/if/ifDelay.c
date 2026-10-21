/**CFile****************************************************************

  FileName    [ifDelay.c]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [FPGA mapping based on priority cuts.]

  Synopsis    [Delay balancing for cut functions.]

  Author      [Alan Mishchenko]
  
  Affiliation [UC Berkeley]

  Date        [Ver. 1.0. Started - November 21, 2006.]

  Revision    [$Id: ifDelay.c,v 1.00 2006/11/21 00:00:00 alanmi Exp $]

***********************************************************************/

#include "if.h"
#include "ifCount.h"
#include "bool/kit/kit.h"

ABC_NAMESPACE_IMPL_START

////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

#define IF_MAX_CUBES 70

////////////////////////////////////////////////////////////////////////
///                     FUNCTION DEFINITIONS                         ///
////////////////////////////////////////////////////////////////////////
 
/**Function*************************************************************

  Synopsis    [Computes the SOP delay using balanced AND decomposition.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
static inline int If_CutMaxCubeSize( Vec_Int_t * vCover, int nVars )
{
    int i, k, Entry, Literal, Count, CountMax = 0;
    Vec_IntForEachEntry( vCover, Entry, i )
    {
        Count = 0;
        for ( k = 0; k < nVars; k++ )
        {
            Literal = (3 & (Entry >> (k << 1)));
            if ( Literal == 1 || Literal == 2 )
                Count++;
        }
        CountMax = Abc_MaxInt( CountMax, Count );
    }
    return CountMax;
}
int If_CutDelaySop( If_Man_t * p, If_Cut_t * pCut )
{
    char * pPerm = If_CutPerm( pCut );
    // delay is calculated using 1+log2(NumFanins)
    static double GateDelays[20] = { 1.00, 1.00, 2.00, 2.58, 3.00, 3.32, 3.58, 3.81, 4.00, 4.17, 4.32, 4.46, 4.58, 4.70, 4.81, 4.91, 5.00, 5.09, 5.17, 5.25 };
    Vec_Int_t * vCover;
    If_Obj_t * pLeaf;
    int i, nLitMax, Delay, DelayMax;
    // mark cut as a user cut
    pCut->fUser = 1;
    if ( pCut->nLeaves == 0 )
        return 0;
    if ( pCut->nLeaves == 1 )
        return (int)If_ObjCutBest(If_CutLeaf(p, pCut, 0))->Delay;
    vCover = Vec_WecEntry( p->vTtIsops[pCut->nLeaves], Abc_Lit2Var(If_CutTruthLit(pCut)) );
    if ( Vec_IntSize(vCover) == 0 )
        return -1;
    // mark the output as complemented
//    vAnds = If_CutDelaySopAnds( p, pCut, vCover, RetValue ^ pCut->fCompl );
    if ( Vec_IntSize(vCover) > p->pPars->nGateSize )
        return -1;
    // set the area cost
    assert( If_CutLeaveNum(pCut) >= 0 && If_CutLeaveNum(pCut) <= 16 );
    // compute the gate delay
    nLitMax = If_CutMaxCubeSize( vCover, If_CutLeaveNum(pCut) );
    if ( Vec_IntSize(vCover) < 2 )
    {
        pCut->Cost = Vec_IntSize(vCover);
        Delay = (int)(GateDelays[If_CutLeaveNum(pCut)] + 0.5);
        DelayMax = 0;
        If_CutForEachLeaf( p, pCut, pLeaf, i )
            DelayMax = Abc_MaxInt( DelayMax, If_ObjCutBest(pLeaf)->Delay + (pPerm[i] = (char)Delay) );
    }
    else
    {
        pCut->Cost = Vec_IntSize(vCover) + 1;
        Delay = (int)(GateDelays[If_CutLeaveNum(pCut)] + GateDelays[nLitMax] + 0.5);
        DelayMax = 0;
        If_CutForEachLeaf( p, pCut, pLeaf, i )
            DelayMax = Abc_MaxInt( DelayMax, If_ObjCutBest(pLeaf)->Delay + (pPerm[i] = (char)Delay) );
    }
    return DelayMax;
}


/**Function*************************************************************

  Synopsis    [Compute pin delays.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int If_CutSopBalancePinDelaysInt( Vec_Int_t * vCover, int * pTimes, word * pFaninRes, int nSuppAll, word * pRes )
{
    word pPinDelsAnd[IF_MAX_FUNC_LUTSIZE], pPinDelsOr[IF_MAX_CUBES];
    int nCounterAnd, pCounterAnd[IF_MAX_FUNC_LUTSIZE];
    int nCounterOr,  pCounterOr[IF_MAX_CUBES];
    int i, k, Entry, Literal, Delay = 0;
    word ResAnd;
    if ( Vec_IntSize(vCover) > IF_MAX_CUBES )
        return -1;
    nCounterOr = 0;
    Vec_IntForEachEntry( vCover, Entry, i )
    { 
        nCounterAnd = 0;
        for ( k = 0; k < nSuppAll; k++ )
        {
            Literal = 3 & (Entry >> (k << 1));
            if ( Literal == 1 || Literal == 2 ) // neg or pos literal
                Delay = If_LogCounterPinDelays( pCounterAnd, &nCounterAnd, pPinDelsAnd, pTimes[k], pFaninRes[k], nSuppAll, 0 );
            else if ( Literal != 0 ) 
                assert( 0 );
        }
        assert( nCounterAnd > 0 );
        ResAnd = If_LogPinDelaysMulti( pPinDelsAnd, nCounterAnd, nSuppAll, 0 );
        Delay = If_LogCounterPinDelays( pCounterOr, &nCounterOr, pPinDelsOr, Delay, ResAnd, nSuppAll, 0 );
    }
    assert( nCounterOr > 0 );
    *pRes = If_LogPinDelaysMulti( pPinDelsOr, nCounterOr, nSuppAll, 0 );
    return Delay;
}
int If_CutSopBalancePinDelaysIntInt( Vec_Int_t * vCover, int * pTimes, int nSuppAll, char * pPerm )
{
    int i, Delay;
    word Res, FaninRes[IF_MAX_FUNC_LUTSIZE];
    for ( i = 0; i < nSuppAll; i++ )
        FaninRes[i] = If_CutPinDelayInit(i);
    Delay = If_CutSopBalancePinDelaysInt( vCover, pTimes, FaninRes, nSuppAll, &Res );
    If_CutPinDelayTranslate( Res, nSuppAll, pPerm );
    return Delay;
}
int If_CutSopBalancePinDelays( If_Man_t * p, If_Cut_t * pCut, char * pPerm )
{
    if ( pCut->nLeaves == 0 ) // const
        return 0;
    if ( pCut->nLeaves == 1 ) // variable
    {
        pPerm[0] = 0;
        return (int)If_ObjCutBest(If_CutLeaf(p, pCut, 0))->Delay;
    }
    else
    {
        Vec_Int_t * vCover;
        int i, pTimes[IF_MAX_FUNC_LUTSIZE];
        vCover = Vec_WecEntry( p->vTtIsops[pCut->nLeaves], Abc_Lit2Var(If_CutTruthLit(pCut)) );
        if ( Vec_IntSize(vCover) == 0 )
            return -1;
        for ( i = 0; i < If_CutLeaveNum(pCut); i++ )
            pTimes[i] = (int)If_ObjCutBest(If_CutLeaf(p, pCut, i))->Delay; 
        return If_CutSopBalancePinDelaysIntInt( vCover, pTimes, If_CutLeaveNum(pCut), pPerm );
    }
}

/**Function*************************************************************

  Synopsis    [Evaluate delay using SOP balancing.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int If_CutSopBalanceEvalInt( Vec_Int_t * vCover, int * pTimes, int * pFaninLits, Vec_Int_t * vAig, int * piRes, int nSuppAll, int * pArea )
{
    int nCounterAnd, pCounterAnd[IF_MAX_FUNC_LUTSIZE], pFaninLitsAnd[IF_MAX_FUNC_LUTSIZE];
    int nCounterOr,  pCounterOr[IF_MAX_CUBES],  pFaninLitsOr[IF_MAX_CUBES];
    int i, k, Entry, Literal, nLits, Delay = 0, iRes = 0;
    if ( Vec_IntSize(vCover) > IF_MAX_CUBES )
        return -1;
    nCounterOr = 0;
    Vec_IntForEachEntry( vCover, Entry, i )
    { 
        nCounterAnd = nLits = 0;
        for ( k = 0; k < nSuppAll; k++ )
        {
            Literal = 3 & (Entry >> (k << 1));
            if ( Literal == 1 ) // neg literal
                nLits++, Delay = If_LogCounterAddAig( pCounterAnd, &nCounterAnd, pFaninLitsAnd, pTimes[k], vAig ? Abc_LitNot(pFaninLits[k]) : -1, vAig, nSuppAll, 0, 0 );
            else if ( Literal == 2 ) // pos literal
                nLits++, Delay = If_LogCounterAddAig( pCounterAnd, &nCounterAnd, pFaninLitsAnd, pTimes[k], vAig ? pFaninLits[k] : -1, vAig, nSuppAll, 0, 0 );
            else if ( Literal != 0 ) 
                assert( 0 );
        }
        assert( nCounterAnd > 0 );
        assert( nLits > 0 );
        if ( vAig )
            iRes = If_LogCreateAndXorMulti( vAig, pFaninLitsAnd, nCounterAnd, nSuppAll, 0 );
        else
            *pArea += nLits == 1 ? 0 : nLits - 1;
        Delay = If_LogCounterAddAig( pCounterOr, &nCounterOr, pFaninLitsOr, Delay, vAig ? Abc_LitNot(iRes) : -1, vAig, nSuppAll, 0, 0 );
    }
    assert( nCounterOr > 0 );
    if ( vAig )
    {
        *piRes = Abc_LitNot( If_LogCreateAndXorMulti( vAig, pFaninLitsOr, nCounterOr, nSuppAll, 0 ) );
        if ( ((vCover->nCap >> 16) & 1) )  // hack to remember complemented attribute
            *piRes = Abc_LitNot( *piRes );
    }
    else       
        *pArea += Vec_IntSize(vCover) == 1 ? 0 : Vec_IntSize(vCover) - 1;
    return Delay;
}
int If_CutSopBalanceEvalIntInt( Vec_Int_t * vCover, int nLeaves, int * pTimes, Vec_Int_t * vAig, int fCompl, int * pArea ) 
{
    int pFaninLits[IF_MAX_FUNC_LUTSIZE];
    int iRes = 0, Res, k;
    if ( vAig )
        for ( k = 0; k < nLeaves; k++ )
            pFaninLits[k] = Abc_Var2Lit(k, 0);
    Res = If_CutSopBalanceEvalInt( vCover, pTimes, pFaninLits, vAig, &iRes, nLeaves, pArea );
    if ( Res == -1 )
        return -1;
    assert( vAig == NULL || Abc_Lit2Var(iRes) == nLeaves + Abc_Lit2Var(Vec_IntSize(vAig)) - 1 );
    if ( vAig )
        Vec_IntPush( vAig, Abc_LitIsCompl(iRes) ^ fCompl );
    assert( vAig == NULL || (Vec_IntSize(vAig) & 1) );
    return Res;
}
int If_CutSopBalanceEval( If_Man_t * p, If_Cut_t * pCut, Vec_Int_t * vAig )
{
    pCut->fUser = 1;
    if ( vAig )
        Vec_IntClear( vAig );
    if ( pCut->nLeaves == 0 ) // const
    {
        assert( Abc_Lit2Var(If_CutTruthLit(pCut)) == 0 );
        if ( vAig )
            Vec_IntPush( vAig, Abc_LitIsCompl(If_CutTruthLit(pCut)) );
        pCut->Cost = 0;
        return 0;
    }
    if ( pCut->nLeaves == 1 ) // variable
    {
        assert( Abc_Lit2Var(If_CutTruthLit(pCut)) == 1 );
        if ( vAig )
            Vec_IntPush( vAig, 0 );
        if ( vAig )
            Vec_IntPush( vAig, Abc_LitIsCompl(If_CutTruthLit(pCut)) );
        pCut->Cost = 0;
        return (int)If_ObjCutBest(If_CutLeaf(p, pCut, 0))->Delay;
    }
    else
    {
        int fVerbose = 0;
        Vec_Int_t * vCover = Vec_WecEntry( p->vTtIsops[pCut->nLeaves], Abc_Lit2Var(If_CutTruthLit(pCut)) );
        int Delay, Area = 0;
        int i, pTimes[IF_MAX_FUNC_LUTSIZE];
        if ( vCover == NULL )
            return -1;
        assert( Vec_IntSize(vCover) > 0 );
        for ( i = 0; i < If_CutLeaveNum(pCut); i++ )
            pTimes[i] = (int)If_ObjCutBest(If_CutLeaf(p, pCut, i))->Delay; 
        Delay = If_CutSopBalanceEvalIntInt( vCover, If_CutLeaveNum(pCut), pTimes, vAig, Abc_LitIsCompl(If_CutTruthLit(pCut)) ^ pCut->fCompl, &Area );
        pCut->Cost = Area;
        if ( fVerbose )
        {
            int Max = 0, Two = 0;
            for ( i = 0; i < If_CutLeaveNum(pCut); i++ )
                Max = Abc_MaxInt( Max, pTimes[i] );
            for ( i = 0; i < If_CutLeaveNum(pCut); i++ )
                if ( pTimes[i] != Max )
                    Two = Abc_MaxInt( Two, pTimes[i] );
            if ( Two + 2 < Max && Max + 3 < Delay )
            {
                for ( i = 0; i < If_CutLeaveNum(pCut); i++ )
                    printf( "%3d ", pTimes[i] );
                for ( ; i < p->pPars->nLutSize; i++ )
                    printf( "    " );
                printf( "-> %3d   ", Delay );
                Dau_DsdPrintFromTruth( If_CutTruthW(p, pCut), If_CutLeaveNum(pCut) );
                Kit_TruthIsopPrintCover( vCover, If_CutLeaveNum(pCut), Abc_LitIsCompl(If_CutTruthLit(pCut)) ^ pCut->fCompl );
                {
                    Vec_Int_t vIsop;
                    int pIsop[64];
                    vIsop.nCap = vIsop.nSize = Abc_Tt6Esop( *If_CutTruthW(p, pCut), pCut->nLeaves, pIsop );
                    vIsop.pArray = pIsop;
                    printf( "ESOP (%d -> %d)\n", Vec_IntSize(vCover), vIsop.nSize );
                    Kit_TruthIsopPrintCover( &vIsop, If_CutLeaveNum(pCut), 0 );
                }
                printf( "\n" );
            }
        }
        return Delay;
    }
}

/**Function*************************************************************

  Synopsis    [Evaluate delay using SOP balancing.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int If_CutLutBalancePinDelays( If_Man_t * p, If_Cut_t * pCut, char * pPerm )
{
    if ( pCut->nLeaves == 0 ) // const
        return 0;
    if ( pCut->nLeaves == 1 ) // variable
    {
        pPerm[0] = 0;
        return (int)If_ObjCutBest(If_CutLeaf(p, pCut, 0))->Delay;
    }
    else
    {
        char * pCutPerm = If_CutDsdPerm( p, pCut );
        int LutSize = p->pPars->pLutStruct[0] - '0';
        int i, Delay, DelayMax = -1;
        assert( (If_CutLeaveNum(pCut) > LutSize) == (pCut->uMaskFunc > 0) );
        for ( i = 0; i < If_CutLeaveNum(pCut); i++ )
        {
            if ( If_CutLeaveNum(pCut) > LutSize && ((pCut->uMaskFunc >> (i << 1)) & 1) )
                pPerm[Abc_Lit2Var((int)pCutPerm[i])] = 2;
            else
                pPerm[Abc_Lit2Var((int)pCutPerm[i])] = 1;
        }
        for ( i = 0; i < If_CutLeaveNum(pCut); i++ )
        {
            Delay = (int)If_ObjCutBest(If_CutLeaf(p, pCut, i))->Delay;
            DelayMax = Abc_MaxInt( DelayMax, Delay + (int)pPerm[i] );
        }
        return DelayMax;
    }
}

/**Function*************************************************************

  Synopsis    [Evaluate delay using SOP balancing.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int If_CutLutBalanceEval( If_Man_t * p, If_Cut_t * pCut )
{
    pCut->fUser = 1;
    pCut->Cost = pCut->nLeaves > 1 ? 1 : 0;
    pCut->uMaskFunc = 0;
    if ( pCut->nLeaves == 0 ) // const
    {
        assert( Abc_Lit2Var(If_CutTruthLit(pCut)) == 0 );
        return 0;
    }
    if ( pCut->nLeaves == 1 ) // variable
    {
        assert( Abc_Lit2Var(If_CutTruthLit(pCut)) == 1 );
        return (int)If_ObjCutBest(If_CutLeaf(p, pCut, 0))->Delay;
    }
    else
    {
        char * pCutPerm = If_CutDsdPerm( p, pCut );
        int LutSize = p->pPars->pLutStruct[0] - '0';
        int i, pTimes[IF_MAX_FUNC_LUTSIZE];
        int DelayMax = -1, nLeafMax = 0;
        unsigned uLeafMask = 0;
        for ( i = 0; i < If_CutLeaveNum(pCut); i++ )
        {
            pTimes[i] = (int)If_ObjCutBest(If_CutLeaf(p, pCut, Abc_Lit2Var((int)pCutPerm[i])))->Delay; 
            if ( DelayMax < pTimes[i] )
                DelayMax = pTimes[i], nLeafMax = 1, uLeafMask = (1 << (i << 1));
            else if ( DelayMax == pTimes[i] )
                nLeafMax++, uLeafMask |= (1 << (i << 1));
        }
        if ( If_CutLeaveNum(pCut) <= LutSize )
            return DelayMax + 1;
        pCut->Cost = 2;
        if ( nLeafMax <= LutSize - 1 )
        {
            pCut->uMaskFunc = If_DsdManCheckXY( p->pIfDsdMan, If_CutDsdLit(p, pCut), LutSize, 1, uLeafMask, 0, 0 );
            if ( pCut->uMaskFunc > 0 )
                return DelayMax + 1;
        }
        pCut->uMaskFunc = If_DsdManCheckXY( p->pIfDsdMan, If_CutDsdLit(p, pCut), LutSize, 1, 0, 0, 0 );
        if ( pCut->uMaskFunc == 0 )
            return -1;
        return DelayMax + 2;
    }
}

/**Function*************************************************************

  Synopsis    [Care-set extraction for LUT decomposition with don't cares.]

  Description [Simulates a window of logic in front of the cut. The window
  inputs start as the cut leaves; the input whose fan-ins add the fewest new
  inputs is replaced by its fan-ins until pPars->nDcWindowSize inputs are
  reached (combinational inputs never expand). Every input pattern of the
  window produces one pattern at the cut leaves, which is marked as a care.
  Marks are cleared through the list of touched objects and the simulation
  buffer is reused, so the cost of a cut depends on the window only.
  Returns 1 if the care set has at least one don't care, 0 if it is all ones.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
#define IF_DC_MAX_WINDOW 16

static inline void If_DcMark( If_Man_t * p, If_Obj_t * pObj )
{
    pObj->fVisit = 1;
    Vec_IntPush( p->vDcTouched, pObj->Id );
}
static inline int If_DcExpandCost( If_Obj_t * pObj )
{
    If_Obj_t * pFanin0 = If_ObjFanin0( pObj ), * pFanin1 = If_ObjFanin1( pObj );
    return ( !pFanin0->fVisit && !If_ObjIsConst1( pFanin0 ) ) + ( !pFanin1->fVisit && !If_ObjIsConst1( pFanin1 ) );
}
// collects the window inputs; the fVisit marks of the touched objects stay set
static int If_DcCreateWindow( If_Man_t * p, If_Cut_t * pCut, Vec_Int_t * vInputs )
{
    If_Obj_t * pLeaf, * pFanin, * pBest;
    int i, Entry, Cost, bestCost, bestPos;
    int nMax = Abc_MinInt( Abc_MaxInt( p->pPars->nDcWindowSize, (int)pCut->nLeaves ), IF_DC_MAX_WINDOW );
    Vec_IntClear( vInputs );
    If_CutForEachLeaf( p, pCut, pLeaf, i )
    {
        If_DcMark( p, pLeaf );
        Vec_IntPush( vInputs, pLeaf->Id );
    }
    while ( 1 )
    {
        bestCost = 100; bestPos = -1; pBest = NULL;
        Vec_IntForEachEntry( vInputs, Entry, i )
        {
            pLeaf = If_ManObj( p, Entry );
            if ( !If_ObjIsAnd( pLeaf ) )
                continue;
            Cost = If_DcExpandCost( pLeaf );
            if ( Cost < bestCost || ( Cost == bestCost && pLeaf->Level > pBest->Level ) )
            {
                bestCost = Cost; bestPos = i; pBest = pLeaf;
            }
            if ( bestCost == 0 )
                break;
        }
        if ( pBest == NULL || Vec_IntSize( vInputs ) - 1 + bestCost > nMax )
            break;
        Vec_IntWriteEntry( vInputs, bestPos, Vec_IntEntryLast( vInputs ) );
        Vec_IntPop( vInputs );
        pFanin = If_ObjFanin0( pBest );
        if ( !pFanin->fVisit && !If_ObjIsConst1( pFanin ) )
        {
            If_DcMark( p, pFanin );
            Vec_IntPush( vInputs, pFanin->Id );
        }
        pFanin = If_ObjFanin1( pBest );
        if ( !pFanin->fVisit && !If_ObjIsConst1( pFanin ) )
        {
            If_DcMark( p, pFanin );
            Vec_IntPush( vInputs, pFanin->Id );
        }
    }
    return Vec_IntSize( vInputs );
}
// collects the window nodes between the inputs and the cut leaves in topological order
static void If_DcCollectNodes_rec( If_Man_t * p, If_Obj_t * pObj, Vec_Int_t * vNodes, int nInputs )
{
    if ( Vec_IntEntry( p->vDcSlot, pObj->Id ) >= 0 )
        return;
    assert( If_ObjIsAnd( pObj ) );
    if ( !If_ObjIsConst1( If_ObjFanin0( pObj ) ) )
        If_DcCollectNodes_rec( p, If_ObjFanin0( pObj ), vNodes, nInputs );
    if ( !If_ObjIsConst1( If_ObjFanin1( pObj ) ) )
        If_DcCollectNodes_rec( p, If_ObjFanin1( pObj ), vNodes, nInputs );
    Vec_IntWriteEntry( p->vDcSlot, pObj->Id, nInputs + Vec_IntSize( vNodes ) );
    Vec_IntPush( vNodes, pObj->Id );
}
int If_ExtractDc( If_Man_t * p, If_Cut_t * pCut, word * pCareSet )
{
    Vec_Int_t * vInputs = p->vDcInputs, * vNodes = p->vDcNodes;
    If_Obj_t * pObj, * pLeaf;
    word * pSim, * pRes, * pF0, * pF1, a, b;
    word * pLeafSim[IF_DC_MAX_WINDOW], LeafWord[IF_DC_MAX_WINDOW];
    int i, k, w, bit, Entry, nInputs, nWords, nLeaves, nLeafWords, nPatBits, nOnes;
    unsigned entry;
    if ( Vec_IntSize( p->vDcSlot ) < If_ManObjNum( p ) )
        Vec_IntFill( p->vDcSlot, If_ManObjNum( p ), -1 );
    // window inputs
    nInputs = If_DcCreateWindow( p, pCut, vInputs );
    Vec_IntForEachEntry( p->vDcTouched, Entry, i )
        If_ManObj( p, Entry )->fVisit = 0;
    Vec_IntClear( p->vDcTouched );
    // window nodes
    Vec_IntForEachEntry( vInputs, Entry, i )
        Vec_IntWriteEntry( p->vDcSlot, Entry, i );
    Vec_IntClear( vNodes );
    If_CutForEachLeaf( p, pCut, pLeaf, i )
        If_DcCollectNodes_rec( p, pLeaf, vNodes, nInputs );
    // window that could not expand (all leaves are combinational inputs): every leaf pattern is reachable
    if ( Vec_IntSize( vNodes ) == 0 )
    {
        Vec_IntForEachEntry( vInputs, Entry, i )
            Vec_IntWriteEntry( p->vDcSlot, Entry, -1 );
        nLeafWords = pCut->nLeaves <= 6 ? 1 : ( 1 << ( pCut->nLeaves - 6 ) );
        Abc_TtFill( pCareSet, nLeafWords );
        return 0;
    }
    // simulation over the window inputs
    nWords = nInputs <= 6 ? 1 : ( 1 << ( nInputs - 6 ) );
    Vec_WrdFill( p->vDcSim, nWords * ( nInputs + Vec_IntSize( vNodes ) ), 0 );
    pSim = Vec_WrdArray( p->vDcSim );
    for ( i = 0; i < nInputs; i++ )
        for ( w = 0; w < nWords; w++ )
            pSim[i * nWords + w] = i < 6 ? s_Truths6[i] : ( ( ( w >> ( i - 6 ) ) & 1 ) ? ~(word)0 : (word)0 );
    Vec_IntForEachEntry( vNodes, Entry, i )
    {
        pObj = If_ManObj( p, Entry );
        pRes = pSim + ( nInputs + i ) * nWords;
        pF0 = If_ObjIsConst1( If_ObjFanin0( pObj ) ) ? NULL : pSim + Vec_IntEntry( p->vDcSlot, If_ObjFanin0( pObj )->Id ) * nWords;
        pF1 = If_ObjIsConst1( If_ObjFanin1( pObj ) ) ? NULL : pSim + Vec_IntEntry( p->vDcSlot, If_ObjFanin1( pObj )->Id ) * nWords;
        for ( w = 0; w < nWords; w++ )
        {
            a = pF0 ? pF0[w] : ~(word)0;
            b = pF1 ? pF1[w] : ~(word)0;
            if ( If_ObjFaninC0( pObj ) ) a = ~a;
            if ( If_ObjFaninC1( pObj ) ) b = ~b;
            pRes[w] = a & b;
        }
    }
    // care set: every window input pattern marks the leaf pattern it produces (word by word, the
    // leaf words of 64 patterns are loaded once and transposed bit by bit)
    nLeaves = (int)pCut->nLeaves;
    nLeafWords = nLeaves <= 6 ? 1 : ( 1 << ( nLeaves - 6 ) );
    memset( pCareSet, 0, sizeof(word) * nLeafWords );
    for ( k = 0; k < nLeaves; k++ )
        pLeafSim[k] = pSim + Vec_IntEntry( p->vDcSlot, pCut->pLeaves[k] ) * nWords;
    nPatBits = nInputs >= 6 ? 64 : ( 1 << nInputs );
    for ( w = 0; w < nWords; w++ )
    {
        for ( k = 0; k < nLeaves; k++ )
            LeafWord[k] = pLeafSim[k][w];
        for ( bit = 0; bit < nPatBits; bit++ )
        {
            entry = 0;
            for ( k = 0; k < nLeaves; k++ )
                entry |= (unsigned)( ( LeafWord[k] >> bit ) & 1 ) << k;
            pCareSet[entry >> 6] |= (word)1 << ( entry & 63 );
        }
    }
    // reset the slots
    Vec_IntForEachEntry( vInputs, Entry, i )
        Vec_IntWriteEntry( p->vDcSlot, Entry, -1 );
    Vec_IntForEachEntry( vNodes, Entry, i )
        Vec_IntWriteEntry( p->vDcSlot, Entry, -1 );
    nOnes = Abc_TtCountOnesVec( pCareSet, nLeafWords );
    return nOnes != ( 1 << pCut->nLeaves );
}

int If_LutDecEval( If_Man_t * p, If_Cut_t * pCut, If_Obj_t * pObj, int optDelay, int fFirst )
{
    pCut->fUser = 1;
    pCut->Cost = pCut->nLeaves > 1 ? 1 : 0;
    pCut->decDelay = 0;
    if ( pCut->nLeaves == 0 ) // const
    {
        assert( Abc_Lit2Var(If_CutTruthLit(pCut)) == 0 );
        return 0;
    }
    if ( pCut->nLeaves == 1 ) // variable
    {
        assert( Abc_Lit2Var(If_CutTruthLit(pCut)) == 1 );
        return (int)If_ObjCutBest(If_CutLeaf(p, pCut, 0))->Delay;
    }

    int LutSize = p->pPars->nLutDecSize;
    int i, leaf_delay;
    int DelayMax = -1, nLeafMax = 0;
    unsigned uLeafMask = 0;
    for ( i = 0; i < If_CutLeaveNum(pCut); i++ )
    {
        leaf_delay = If_ObjCutBest(If_CutLeaf(p, pCut, i))->Delay;

        if ( DelayMax < leaf_delay )
        {
            DelayMax = leaf_delay;
            nLeafMax = 1;
            uLeafMask = (1 << i);
        }
        else if ( DelayMax == leaf_delay )
        {
            nLeafMax++;
            uLeafMask |= (1 << i);
        }
    }
    if ( If_CutLeaveNum(pCut) <= LutSize )
    {
        pCut->decDelay = ( 1 << LutSize ) - 1;
        return DelayMax + 1;
    }

    /* compute the decomposition */
    int use_late_arrival = 0;
    unsigned cost = 1;

    if ( !fFirst )
    {
        if ( optDelay )
        {
            /* checks based on delay: must be better than the previous best cut */
            use_late_arrival = DelayMax + 2 >= If_ObjCutBest(pObj)->Delay;
        }
        else
        {
            /* checks based on delay: look at the required time */
            use_late_arrival = DelayMax + 2 > pObj->Required + p->fEpsilon;
        }
    }

    /* Too many late-arriving signals */
    if ( nLeafMax == LutSize )
    {
        if ( use_late_arrival )
        {
            /* unfeasible decomposition */
            pCut->Cost = IF_COST_MAX;
            return ABC_INFINITY;
        }
        else
        {
            /* remove critical signals as not needed */
            uLeafMask = 0;
        }
    }

    /* returns the delay of the decomposition */
    word *pTruth = If_CutTruthW( p, pCut );
    int val;
    if ( p->pPars->fUserLutDecDc )
    {
        // evaluate without don't cares (the baseline of -Z) and, if the care set has a don't
        // care, with them; the decomposition with don't cares is taken only when it is strictly
        // better (feasible where the plain one is not, lower delay, or fewer LUTs) and only then
        // is the care set stored for Abc_DecRecordToHop
        word pCareSet[1 << ( IF_DC_MAX_WINDOW - 6 )];
        unsigned uLeafMaskDc = uLeafMask, costDc = 1;
        word Key[9];
        int * pKey = (int *)Key;
        int valDc, iLeafSet, iCare;
        abctime clk;
        pCut->fDcHasDc = pCut->fDcOnly = pCut->fDcCheaper = 0;
        pCut->iCutCare = -1;
        // the column multiplicity with don't cares is a heuristic whose result depends on the
        // variable order, so a later run with the free set as profile may not find the
        // decomposition again (mem_ctrl, K = 10, 2026-10-20); Abc_DecRecordToHop therefore replays
        // the evaluation with the same inputs
        pCut->decLate = uLeafMask;
        pCut->fDcTryNoLate = !use_late_arrival;
        p->nDcCutsEval++;
        clk = Abc_Clock();
        // the care set depends on the leaves only: look the leaf set up before extracting
        memset( Key, 0, sizeof(Key) );
        pKey[0] = pCut->nLeaves;
        for ( i = 0; i < (int)pCut->nLeaves; i++ )
            pKey[1 + i] = pCut->pLeaves[i];
        iLeafSet = Vec_MemHashInsert( p->vDcLeafMem, Key );
        if ( iLeafSet < Vec_IntSize( p->vDcLeafCare ) )
        {
            iCare = Vec_IntEntry( p->vDcLeafCare, iLeafSet );
            pCut->fDcHasDc = iCare >= 0;
            if ( pCut->fDcHasDc )
                Abc_TtCopy( pCareSet, Vec_MemReadEntry( p->vTtMem[pCut->nLeaves], Abc_Lit2Var( iCare ) ), p->nTruth6Words[pCut->nLeaves], 0 );
        }
        else
        {
            assert( iLeafSet == Vec_IntSize( p->vDcLeafCare ) );
            pCut->fDcHasDc = If_ExtractDc( p, pCut, pCareSet );
            iCare = pCut->fDcHasDc ? Abc_Var2Lit( Vec_MemHashInsert( p->vTtMem[pCut->nLeaves], pCareSet ), 0 ) : -1;
            Vec_IntPush( p->vDcLeafCare, iCare );
        }
        p->timeDcExtract += Abc_Clock() - clk;
        clk = Abc_Clock();
        val = acd_evaluate( pTruth, pCut->nLeaves, LutSize, &uLeafMask, &cost, !use_late_arrival );
        p->timeDcPlain += Abc_Clock() - clk;
        if ( pCut->fDcHasDc )
        {
            p->nDcCutsWithDc++;
            clk = Abc_Clock();
            valDc = acd_dc_evaluate( pTruth, pCareSet, pCut->nLeaves, LutSize, &uLeafMaskDc, &costDc, !use_late_arrival );
            p->timeDcEval += Abc_Clock() - clk;
            if ( valDc >= 0 && ( val < 0 || valDc < val || ( valDc == val && costDc < cost ) ) )
            {
                if ( val < 0 )
                {
                    p->nDcCutsOnly++;
                    pCut->fDcOnly = 1;
                }
                else
                {
                    p->nDcCutsCheaper++;
                    pCut->fDcCheaper = 1;
                }
                val = valDc;
                cost = costDc;
                uLeafMask = uLeafMaskDc;
                pCut->iCutCare = iCare;
            }
        }
    }
    else
    {
        val = acd_evaluate( pTruth, pCut->nLeaves, LutSize, &uLeafMask, &cost, !use_late_arrival );
    }

    /* not feasible decomposition */
    pCut->decDelay = uLeafMask;
    if ( val < 0 )
    {
        pCut->Cost = IF_COST_MAX;
        return ABC_INFINITY;
    }

    pCut->Cost = cost;

    return DelayMax + val;
}

int If_Lut2DecEval( If_Man_t * p, If_Cut_t * pCut, If_Obj_t * pObj, int optDelay, int fFirst )
{
    pCut->fUser = 1;
    pCut->Cost = pCut->nLeaves > 1 ? 1 : 0;
    pCut->decDelay = 0;
    if ( pCut->nLeaves == 0 ) // const
    {
        assert( Abc_Lit2Var(If_CutTruthLit(pCut)) == 0 );
        return 0;
    }
    if ( pCut->nLeaves == 1 ) // variable
    {
        assert( Abc_Lit2Var(If_CutTruthLit(pCut)) == 1 );
        return (int)If_ObjCutBest(If_CutLeaf(p, pCut, 0))->Delay;
    }

    int LutSize = p->pPars->nLutDecSize;
    int i, leaf_delay;
    int DelayMax = -1, nLeafMax = 0;
    unsigned uLeafMask = 0;
    for ( i = 0; i < If_CutLeaveNum(pCut); i++ )
    {
        leaf_delay = If_ObjCutBest(If_CutLeaf(p, pCut, i))->Delay;

        if ( DelayMax < leaf_delay )
        {
            DelayMax = leaf_delay;
            nLeafMax = 1;
            uLeafMask = (1 << i);
        }
        else if ( DelayMax == leaf_delay )
        {
            nLeafMax++;
            uLeafMask |= (1 << i);
        }
    }
    if ( If_CutLeaveNum(pCut) <= LutSize )
    {
        pCut->decDelay = ( 1 << LutSize ) - 1;
        return DelayMax + 1;
    }

    /* compute the decomposition */
    int use_late_arrival = 0;
    unsigned cost = 1;

    if ( !fFirst )
    {
        if ( optDelay )
        {
            /* checks based on delay: must be better than the previous best cut */
            use_late_arrival = DelayMax + 2 >= If_ObjCutBest(pObj)->Delay;
        }
        else
        {
            /* checks based on delay: look at the required time */
            use_late_arrival = DelayMax + 2 > pObj->Required + p->fEpsilon;
        }
    }

    /* Too many late-arriving signals */
    if ( nLeafMax == LutSize && use_late_arrival )
    {
        /* unfeasible decomposition */
        pCut->Cost = IF_COST_MAX;
        return ABC_INFINITY;
    }

    if ( !use_late_arrival )
    {
        uLeafMask = 0;
    }

    /* returns the delay of the decomposition */
    word *pTruth = If_CutTruthW( p, pCut );
    int val = acd2_evaluate( pTruth, pCut->nLeaves, LutSize, &uLeafMask, &cost, !use_late_arrival );

    /* not feasible decomposition */
    pCut->decDelay = uLeafMask;
    if ( val < 0 )
    {
        pCut->Cost = IF_COST_MAX;
        return ABC_INFINITY;
    }

    pCut->Cost = 2;
    return DelayMax + val;
}

int If_LutDecReEval( If_Man_t * p, If_Cut_t * pCut )
{
    // pCut->fUser = 1;

    if ( pCut->nLeaves == 0 ) // const
    {
        assert( Abc_Lit2Var(If_CutTruthLit(pCut)) == 0 );
        return 0;
    }
    if ( pCut->nLeaves == 1 ) // variable
    {
        assert( Abc_Lit2Var(If_CutTruthLit(pCut)) == 1 );
        return (int)If_ObjCutBest(If_CutLeaf(p, pCut, 0))->Delay;
    }

    // int LutSize = p->pPars->pLutStruct[0] - '0';
    int i, leaf_delay;
    int DelayMax = -1;
    for ( i = 0; i < If_CutLeaveNum(pCut); i++ )
    {
        leaf_delay = If_ObjCutBest(If_CutLeaf(p, pCut, i))->Delay;
        leaf_delay += ( ( pCut->decDelay >> i ) & 1 ) == 0 ? 2 : 1;
        DelayMax = Abc_MaxInt( leaf_delay, DelayMax );
    }

    return DelayMax;
}

float If_LutDecPinRequired( If_Man_t * p, If_Cut_t * pCut, int i, float required )
{
    if ( pCut->nLeaves == 0 ) // const
    {
        assert( Abc_Lit2Var(If_CutTruthLit(pCut)) == 0 );
        return required;
    }
    if ( pCut->nLeaves == 1 ) // variable
    {
        assert( Abc_Lit2Var(If_CutTruthLit(pCut)) == 1 );
        return 0;
    }

    return ( ( pCut->decDelay >> i ) & 1 ) == 0 ? 2 : 1;
}

/*
int If_CutLutBalanceEval( If_Man_t * p, If_Cut_t * pCut )
{
    char pPerm[16];
    int Delay2, Delay = If_CutLutBalanceEvalInt( p, pCut );
    if ( Delay == -1 )
        return -1;
    Delay2 = If_CutLutBalancePinDelays( p, pCut, pPerm );
    if ( Delay2 != Delay )
    {
        int s = 0;
        char * pCutPerm = If_CutDsdPerm( p, pCut );
        If_DsdManPrintNode( p->pIfDsdMan, If_CutDsdLit(p, pCut) );        Dau_DecPrintSet( pCut->uMaskFunc, pCut->nLeaves, 1 );
        Kit_DsdPrintFromTruth( If_CutTruthUR(p, pCut), pCut->nLeaves ); printf( "\n" );
        for ( s = 0; s < pCut->nLeaves; s++ )
//            printf( "%d ", (int)If_ObjCutBest(If_CutLeaf(p, pCut, Abc_Lit2Var((int)pCutPerm[s])))->Delay );
            printf( "%d ", (int)If_ObjCutBest(If_CutLeaf(p, pCut, s))->Delay );
        printf( "\n" );

        Delay  = If_CutLutBalanceEvalInt( p, pCut );
        Delay2 = If_CutLutBalancePinDelays( p, pCut, pPerm );
    }

    return Delay;
}
*/

////////////////////////////////////////////////////////////////////////
///                       END OF FILE                                ///
////////////////////////////////////////////////////////////////////////


ABC_NAMESPACE_IMPL_END

