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

#include <stdint.h>
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

static inline int If_NodeGetLeafCostOne( If_Obj_t * pObj )
{
    int Cost;

    assert( pObj->fVisit == 1 );  // must be in the current cone

    if ( If_ObjIsCi( pObj ) )
        return 999;

    Cost = (!If_ObjFanin0(pObj)->fVisit) + (!If_ObjFanin1(pObj)->fVisit);

    if (pObj->nFanouts > 100000)
    {
        return 999;
    }

    return Cost;
}

int If_ManCreateWindowOld( If_Man_t * p, If_Obj_t * pObj, If_Cut_t * pCut, Vec_Int_t * vNodes, int maxLeaves )
{
    If_Obj_t * pLeaf, * pFanin, * pBest = NULL;
    int i, bestCost, bestPos, cost;

    // Reset fVisit across the network
    If_ManCleanMarkV( p );

    // Initialize window with fanins of pObj
    p->pWindow->nLeaves = 0;
    Vec_IntPush( vNodes, pObj->Id );

    pFanin = If_ObjFanin0( pObj );
    pFanin->fVisit = 1;
    p->pWindow->pLeaves[ p->pWindow->nLeaves++ ] = pFanin->Id;
    Vec_IntPush( vNodes, pFanin->Id );

    pFanin = If_ObjFanin1( pObj );
    pFanin->fVisit = 1;
    p->pWindow->pLeaves[ p->pWindow->nLeaves++ ] = pFanin->Id;
    Vec_IntPush( vNodes, pFanin->Id );

    // Greedily expand the window
    while ( 1 )
    {
        bestCost = 100;
        bestPos = -1;
        pBest = NULL;

        for ( i = 0; i < p->pWindow->nLeaves; ++i )
        {
            pLeaf = If_ManObj( p, p->pWindow->pLeaves[i] );
            cost = If_NodeGetLeafCostOne( pLeaf );

            if ( cost < bestCost || (cost == bestCost && pBest && pLeaf->Level > pBest->Level) )
            {
                bestCost = cost;
                bestPos = i;
                pBest = pLeaf;
            }

            if ( bestCost == 0 )
                break;
        }

        if ( pBest == NULL || p->pWindow->nLeaves - 1 + bestCost > maxLeaves )
            break;

        // Remove pBest from leaves
        for ( i = bestPos; i < p->pWindow->nLeaves - 1; ++i )
            p->pWindow->pLeaves[i] = p->pWindow->pLeaves[i + 1];
        p->pWindow->nLeaves--;

        // Add fanins of pBest
        pFanin = If_ObjFanin0( pBest );
        if ( !pFanin->fVisit )
        {
            pFanin->fVisit = 1;
            p->pWindow->pLeaves[ p->pWindow->nLeaves++ ] = pFanin->Id;
            Vec_IntPush( vNodes, pFanin->Id );
        }

        pFanin = If_ObjFanin1( pBest );
        if ( !pFanin->fVisit )
        {
            pFanin->fVisit = 1;
            p->pWindow->pLeaves[ p->pWindow->nLeaves++ ] = pFanin->Id;
            Vec_IntPush( vNodes, pFanin->Id );
        }

        assert( p->pWindow->nLeaves <= maxLeaves );
    }

    // Check if all leaves in pCut are contained
    for ( i = 0; i < (int)pCut->nLeaves; ++i )
    {
        If_Obj_t * pCutLeaf = If_ManObj( p, pCut->pLeaves[i] );
        if ( !pCutLeaf->fVisit )
        {
            return 0;
        }
    }

    return 1;
}

static const word s_Proj6[6] = {
        0xaaaaaaaaaaaaaaaaULL, // x0 = 101010...
        0xccccccccccccccccULL, // x1 = 11001100...
        0xf0f0f0f0f0f0f0f0ULL, // x2 = 11110000...
        0xff00ff00ff00ff00ULL, // x3 = 8-8 alternating 0s/1s
        0xffff0000ffff0000ULL, // x4 = 16-16 alternating
        0xffffffff00000000ULL  // x5 = 32-32 alternating
};

void If_ManTruthCreateVar( word* pTruth, int nVars, int varIdx, int fCompl )
{
    assert( varIdx < nVars );
    assert( nVars <= 16 );

    const int nBits  = (1 << nVars);
    const word nWords = ( nVars <= 6 ) ? 1 : ( 1 << ( nVars - 6 ) );

    if ( varIdx < 6 && nVars <= 6 )
    {
        pTruth[0] = fCompl ? ~s_Proj6[varIdx] : s_Proj6[varIdx];
        // Mask unused bits in the last word
        if ( nBits < 64 )
        {
            pTruth[0] &= ~(~(word)0 << nBits);
        }
        return;
    }

    if ( varIdx < 6 )
    {
        word pattern = fCompl ? ~s_Proj6[varIdx] : s_Proj6[varIdx];
        for ( int i = 0; i < nWords; ++i )
            pTruth[i] = pattern;
        return;
    }

    // varIdx ≥ 6 → alternating blocks
    const int c = 1 << ( varIdx - 6 );
    const word zero = UINT64_C( 0 );
    const word one = ~zero;
    ABC_UINT64_T block = UINT64_C( 0 );

    while ( block < nWords )
    {
        for ( int i = 0; i < c; ++i )
        {
            pTruth[block++] = fCompl ? one : zero;
        }
        for ( int i = 0; i < c; ++i )
        {
            pTruth[block++] = fCompl ? zero : one;
        }
    }
}

void If_ManSimulateWindowOld( If_Man_t* p, If_Obj_t * pRoot, If_Cut_t* pCut, Vec_Int_t* vNodes, Vec_Ptr_t* vSimTts )
{
    int Entry, i;
    int nWords = ( p->pWindow->nLeaves <= 6 ) ? 1 : ( 1 << ( p->pWindow->nLeaves - 6 ) );

    // Initialize PI values for each leaf of the window
    for ( i = 0; i < p->pWindow->nLeaves; ++i )
    {
        int piId = p->pWindow->pLeaves[i];
        word* pTruth = ABC_ALLOC( word, nWords );
        If_ManTruthCreateVar( pTruth, p->pWindow->nLeaves, i, 0 );
        Vec_PtrWriteEntry( vSimTts, piId, pTruth );
    }

    // Simulate internal nodes in DFS order
    Vec_IntForEachEntry( vNodes, Entry, i )
    {
        If_Obj_t* pObj = If_ManObj( p, Entry );

        word* pFan0Orig = (word*)Vec_PtrEntry( vSimTts, If_ObjFanin0(pObj)->Id );
        word* pFan1Orig = (word*)Vec_PtrEntry( vSimTts, If_ObjFanin1(pObj)->Id );

        // Allocate result for this node
        word* pResult = ABC_ALLOC( word, nWords );

        // Handle complemented fanins using temp copies
        word* pFan0 = ABC_ALLOC( word, nWords );
        word* pFan1 = ABC_ALLOC( word, nWords );
        Abc_TtCopy( pFan0, pFan0Orig, nWords, If_ObjFaninC0(pObj) );
        Abc_TtCopy( pFan1, pFan1Orig, nWords, If_ObjFaninC1(pObj) );

        // Compute AND and store result
        Abc_TtAnd( pResult, pFan0, pFan1, nWords, 0 );
        Vec_PtrWriteEntry( vSimTts, pObj->Id, pResult );

        ABC_FREE( pFan0 );
        ABC_FREE( pFan1 );
    }
}

void If_ManSortWindowNodes( If_Man_t * p, Vec_Int_t * vNodes )
{
    int i, Entry;
    If_Obj_t * pObj;

    If_ManCleanMarkV( p );

    Vec_Int_t * vTemp = Vec_IntAlloc( Vec_IntSize( vNodes ) );

    for ( i = 0; i < (int)p->pWindow->nLeaves; ++i )
    {
        pObj = If_ManObj( p, p->pWindow->pLeaves[i] );
        pObj->fVisit = 1;
        // Vec_IntPush( vTemp, pObj->Id );
    }

    int nRemain = Vec_IntSize( vNodes ) - p->pWindow->nLeaves;
    int nAdded = 1;

    while ( nRemain && nAdded )
    {
        nAdded = 0;
        Vec_IntForEachEntry( vNodes, Entry, i )
        {
            pObj = If_ManObj( p, Entry );
            if ( pObj->fVisit )
                continue;

            If_Obj_t * pFanin0 = If_ObjFanin0( pObj );
            If_Obj_t * pFanin1 = If_ObjFanin1( pObj );

            if ( pFanin0->fVisit && pFanin1->fVisit )
            {
                pObj->fVisit = 1;
                Vec_IntPush( vTemp, Entry );
                ++nAdded;
            }
        }
        nRemain -= nAdded;
    }

    if ( nRemain > 0 )
        printf( "Warning: Not all window nodes were sorted due to unresolved fanins.\n" );

    // Step 4: Overwrite vNodes with sorted result
    Vec_IntClear( vNodes );
    Vec_IntForEachEntry( vTemp, Entry, i )
        Vec_IntPush( vNodes, Entry );

    Vec_IntFree( vTemp );
}

word If_ManGetBit( const word * pTt, int index )
{
    return ( pTt[index >> 6] >> ( index & 0x3f ) ) & 0x1;
}

void If_ManSetBit( word * pTt, int index )
{
    pTt[index >> 6] |= UINT64_C( 1 ) << ( index & 0x3f );
}

void If_ManComputeCareSet( If_Man_t * p, If_Cut_t * pCut, Vec_Ptr_t * vSimTts, word * pCareSet, int nVars )
{
    int i, j;
    If_Obj_t * pLeaf;

    for ( i = 0; i < ( 1u << nVars ); ++i )  // ✅ correct simulation domain
    {
        word entry = 0;

        If_CutForEachLeaf( p, pCut, pLeaf, j )
        {
            const word * NodeTt = (const word *)Vec_PtrEntry( vSimTts, pLeaf->Id );
            entry |= ((word)If_ManGetBit( NodeTt, i )) << j;
        }

        If_ManSetBit( pCareSet, (int)entry );
    }
}

int If_ExtractDcOld( If_Man_t * p, If_Cut_t * pCut, If_Obj_t * pObj, word * pCareSet )
{
    int maxLeaves = 12;
    int ret = 0;

    // Allocate the node vector
    Vec_Int_t * vNodes = Vec_IntAlloc( 64 ); // initial capacity; growable

    // Clear the reusable window
    memset( p->pWindow, 0, sizeof(If_Cut_t) + sizeof(int) * ( maxLeaves + p->nPermWords ) );

    // Create and simulate the window
    if ( If_ManCreateWindowOld( p, pObj, pCut, vNodes, maxLeaves ) )
    {
        // printf("DCs get computed\n");
        If_ManSortWindowNodes( p, vNodes );
        Vec_Ptr_t* vSimTts = Vec_PtrStart( If_ManObjNum(p) );
        If_ManSimulateWindowOld( p, pObj, pCut, vNodes, vSimTts );
        If_ManComputeCareSet( p, pCut, vSimTts, pCareSet, pCut->nLeaves );
        for ( int i = 0; i < If_ManObjNum(p); ++i )
        {
            word* pE = (word*)Vec_PtrEntry( vSimTts, i );
            if ( pE )
                ABC_FREE( pE );
        }
        Vec_PtrFree( vSimTts );
        ret = 1;
    }
    else
    {
        int nWords = ( pCut->nLeaves <= 6 ) ? 1 : ( 1 << ( pCut->nLeaves - 6 ) );
        for ( int i = 0; i < nWords; ++i )
        {
            pCareSet[i] = 0xffffffffffffffffULL;
        }
    }

    // Free node vector
    Vec_IntFree( vNodes );

    return ret;
}

void If_ManSimulateWindow( If_Man_t* p, Vec_Int_t* vInputs, Vec_Int_t* vNodes, Vec_Ptr_t* vSimTts )
{
    int Entry, i;
    int nWords = ( Vec_IntSize(vInputs) <= 6 ) ? 1 : ( 1 << ( Vec_IntSize(vInputs) - 6 ) );

    // Initialize PI values for each leaf of the window
    Vec_IntForEachEntry(vInputs, Entry, i)
    {
        word* pTruth = ABC_ALLOC( word, nWords );
        If_ManTruthCreateVar( pTruth, Vec_IntSize(vInputs), i, 0 );
        Vec_PtrWriteEntry( vSimTts, Entry, pTruth );
    }

    // Simulate internal nodes in topo order
    Vec_IntForEachEntry( vNodes, Entry, i )
    {
        If_Obj_t* pObj = If_ManObj( p, Entry );

        if ( If_ObjIsCi(pObj) )
            printf("This does not make sense");

        word* pFan0Orig = (word*)Vec_PtrEntry( vSimTts, If_ObjFanin0(pObj)->Id );
        word* pFan1Orig = (word*)Vec_PtrEntry( vSimTts, If_ObjFanin1(pObj)->Id );

        // Allocate result for this node
        word* pResult = ABC_ALLOC( word, nWords );

        // Handle complemented fanins using temp copies
        word* pFan0 = ABC_ALLOC( word, nWords );
        word* pFan1 = ABC_ALLOC( word, nWords );
        Abc_TtCopy( pFan0, pFan0Orig, nWords, If_ObjFaninC0(pObj) );
        Abc_TtCopy( pFan1, pFan1Orig, nWords, If_ObjFaninC1(pObj) );

        // Compute AND and store result
        Abc_TtAnd( pResult, pFan0, pFan1, nWords, 0 );
        Vec_PtrWriteEntry( vSimTts, pObj->Id, pResult );

        ABC_FREE( pFan0 );
        ABC_FREE( pFan1 );
    }
}

void If_ManCollectNodesRec( If_Man_t *p, If_Obj_t *pNode, Vec_Int_t *vNodes )
{
    If_Obj_t *pFanin;

    if ( pNode->fVisit == 1 )
        return;

    pNode->fVisit = 1;

    // for each fan-in recursively collect the nodes
    pFanin = If_ObjFanin0( pNode );
    if (If_ObjIsConst1( pFanin ))
        return;
    If_ManCollectNodesRec( p, pFanin, vNodes );
    pFanin = If_ObjFanin1( pNode );
    if (If_ObjIsConst1( pFanin ))
        return;
    If_ManCollectNodesRec( p, pFanin, vNodes );

    Vec_IntPush( vNodes, pNode->Id );
}

void If_ManCollectNodes( If_Man_t *p, If_Cut_t *pCut, Vec_Int_t *vInputs, Vec_Int_t *vNodes )
{
    int i, j, Entry;
    If_Obj_t *pInput, *pLeaf;

    If_ManCleanMarkV( p );

    Vec_IntForEachEntry( vInputs, Entry, i )
    {
        pInput = If_ManObj(p, Entry);
        pInput->fVisit = 1;
    }

    If_CutForEachLeaf( p, pCut, pLeaf, j )
    {
        If_ManCollectNodesRec( p, pLeaf, vNodes );
    }

    If_CutForEachLeaf( p, pCut, pLeaf, j )
    {
        if ( !pLeaf->fVisit )
        {
            pLeaf->fVisit = 1;
            Vec_IntPush( vNodes, pLeaf->Id );
        }
    }
}

int If_ManCreateWindowMin( If_Man_t * p, If_Cut_t * pCut, Vec_Int_t * vInputs )
{
    int i;
    If_Obj_t * pLeaf, * pFanin;

    If_ManCleanMarkV( p );

    If_CutForEachLeaf( p, pCut, pLeaf, i )
    {
        if ( If_ObjIsCi( pLeaf ) )
        {
            Vec_IntPush( vInputs, pLeaf->Id );
            continue;
        }

        pFanin = If_ObjFanin0( pLeaf );
        if ( pFanin && !If_ObjIsConst1( pFanin ) && !pFanin->fVisit )
        {
            pFanin->fVisit = 1;
            Vec_IntPush( vInputs, pFanin->Id );
        }

        pFanin = If_ObjFanin1( pLeaf );
        if ( pFanin && !If_ObjIsConst1( pFanin ) && !pFanin->fVisit )
        {
            pFanin->fVisit = 1;
            Vec_IntPush( vInputs, pFanin->Id );
        }
    }

    return 1;
}

int If_ManCreateWindow( If_Man_t * p, If_Cut_t * pCut, Vec_Int_t * vInputs )
{
    int i, bestCost, bestPos, cost, Entry;
    If_Obj_t * pLeaf, * pFanin, * pBest;

    int maxLeaves = 12;

    If_ManCleanMarkV( p );

    // push all Cut Leaves to vInputs
    If_CutForEachLeaf( p, pCut, pLeaf, i )
    {
        pLeaf->fVisit = 1;
        Vec_IntPush( vInputs, pLeaf->Id );
    }

    // Greedily expand the window
    while ( 1 )
    {
        bestCost = 100;
        bestPos = -1;
        pBest = NULL;

        Vec_IntForEachEntry( vInputs, Entry, i )
        {
            pLeaf = If_ManObj( p, Entry );
            cost = If_NodeGetLeafCostOne( pLeaf );

            if ( cost < bestCost || (cost == bestCost && pBest && pLeaf->Level > pBest->Level) )
            {
                bestCost = cost;
                pBest = pLeaf;
                bestPos = i;
            }

            if ( bestCost == 0 )
                break;
        }

        assert(bestPos >= 0 && bestPos < Vec_IntSize(vInputs) || pBest == NULL);
        if ( pBest == NULL || Vec_IntSize(vInputs) - 1 + bestCost > maxLeaves )
            break;

        // Move to last and delete last
        if ( bestPos != Vec_IntSize(vInputs) - 1 )
            Vec_IntWriteEntry(vInputs, bestPos, Vec_IntEntry(vInputs, Vec_IntSize(vInputs) - 1));
        vInputs->nSize--;

        // Add fanins of pBest to vInputs
        pFanin = If_ObjFanin0( pBest );
        if ( pFanin && pFanin->fVisit != 1 && !If_ObjIsConst1(pFanin) )
        {
            pFanin->fVisit = 1;
            Vec_IntPush( vInputs, pFanin->Id );
        }

        pFanin = If_ObjFanin1( pBest );
        if ( pFanin && pFanin->fVisit != 1 && !If_ObjIsConst1(pFanin) )
        {
            pFanin->fVisit = 1;
            Vec_IntPush( vInputs, pFanin->Id );
        }
    }

    // Check if all leaves in pCut are contained
    for ( i = 0; i < (int)pCut->nLeaves; ++i )
    {
        If_Obj_t * pCutLeaf = If_ManObj( p, pCut->pLeaves[i] );
        if ( !pCutLeaf->fVisit )
        {
            return 0;
        }
    }

    return 1;
}

int If_ExtractDcMin( If_Man_t * p, If_Cut_t * pCut, word * pCareSet )
{
    // Allocate the node vector
    Vec_Int_t * vInputs = Vec_IntAlloc( 64 ); // initial capacity; growable
    Vec_Int_t * vNodes = Vec_IntAlloc( 64 ); // initial capacity; growable

    // Collect the Window Inputs
    If_ManCreateWindow( p, pCut, vInputs );

    // Collect the window nodes in topological order
    If_ManCollectNodes( p, pCut, vInputs, vNodes );

    // Simulate the window
    Vec_Ptr_t* vSimTts = Vec_PtrStart( If_ManObjNum(p) );
    If_ManSimulateWindow( p, vInputs, vNodes, vSimTts );
    If_ManComputeCareSet( p, pCut, vSimTts, pCareSet, Vec_IntSize(vInputs) );

    // Free the data structures used
    for ( int i = 0; i < If_ManObjNum(p); ++i )
    {
        word* pE = (word*)Vec_PtrEntry( vSimTts, i );
        if ( pE )
            ABC_FREE( pE );
    }
    Vec_PtrFree( vSimTts );
    Vec_IntFree( vInputs );
    Vec_IntFree( vNodes );

    /*int nWords = ( pCut->nLeaves <= 6 ) ? 1 : ( 1 << ( pCut->nLeaves - 6 ) );
    int allOnes = 1;

    for ( int i = 0; i < nWords; ++i )
    {
        if ( pCareSet[i] != ~(word)0 )
        {
            allOnes = 0;
            break;
        }
    }

    if ( !allOnes )
        printf( "Dont cares found.\n" );*/

    return 1;
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
        // acd pointer
        unsigned uLeafMaskDc = uLeafMask;
        unsigned costDC = cost;
        // extract the care set
        int nWords = ( pCut->nLeaves <= 6 ) ? 1 : ( 1 << ( pCut->nLeaves - 6 ) );
        word* pCareSet = ABC_ALLOC( word, nWords );
        memset( pCareSet, 0, sizeof(word) * nWords );
        If_ExtractDcMin(p, pCut, pCareSet);

        // save Care set for cut
        int csId = Vec_MemHashInsert(p->vTtMem[pCut->nLeaves], pCareSet);
        pCut->iCutCs = Abc_Var2Lit(csId, 0);

        const int num_blocks = ( pCut->nLeaves <= 6 ) ? 1 : ( 1 << ( pCut->nLeaves - 6 ) );
        int allOnes = 1;
        for ( int i = 0; i < num_blocks; ++i )
        {
            if ( pCareSet[i] != ~(word)0 )
            {
                allOnes = 0;
                break;
            }
        }

        if ( !allOnes )
        {
            /*printf("CareSet: \n");
            for ( int i = 0; i < num_blocks; ++i )
            {
                printf("Block %i: %lu\n", i, pCareSet[i]);
            }*/
            printf( "Dont cares evaluated.\n" );
        }

        val = acd_dc_evaluate( pTruth, pCareSet, pCut->nLeaves, LutSize, &uLeafMask, &cost, !use_late_arrival );
        if ( val != -1 )
        {
            int val2 = acd_evaluate( pTruth, pCut->nLeaves, LutSize, &uLeafMaskDc, &costDC, !use_late_arrival );
            if ( val2 == -1 )
            {
                printf("Decomposition found only using DCs\n");
            }
        }
        ABC_FREE( pCareSet );
    }
    else
    {
        val = acd_evaluate( pTruth, pCut->nLeaves, LutSize, &uLeafMask, &cost, !use_late_arrival );
    }
    //word *pCS= If_CutCsW( p, pCut );

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

