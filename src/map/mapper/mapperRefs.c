/**CFile****************************************************************

  FileName    [mapperRefs.c]

  PackageName [MVSIS 1.3: Multi-valued logic synthesis system.]

  Synopsis    [Generic technology mapping engine.]

  Author      [MVSIS Group]
  
  Affiliation [UC Berkeley]

  Date        [Ver. 2.0. Started - June 1, 2004.]

  Revision    [$Id: mapperRefs.h,v 1.0 2003/09/08 00:00:00 alanmi Exp $]

***********************************************************************/

#include "mapperInt.h"

ABC_NAMESPACE_IMPL_START


////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
///                     FUNCTION DEFINITIONS                         ///
////////////////////////////////////////////////////////////////////////

/**Function*************************************************************

  Synopsis    [Reads the actual reference counter of a phase.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int Map_NodeReadRefPhaseAct( Map_Node_t * pNode, int fPhase )
{
    assert( !Map_IsComplement(pNode) );
    if ( pNode->pCutBest[0] && pNode->pCutBest[1] ) // both assigned
        return pNode->nRefAct[fPhase];
    assert( pNode->pCutBest[0] || pNode->pCutBest[1] ); // at least one assigned
    return pNode->nRefAct[2];
}

/**Function*************************************************************

  Synopsis    [Reads the estimated reference counter of a phase.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
float Map_NodeReadRefPhaseEst( Map_Node_t * pNode, int fPhase )
{
    assert( !Map_IsComplement(pNode) );
    if ( pNode->pCutBest[0] && pNode->pCutBest[1] ) // both assigned
        return pNode->nRefEst[fPhase];
    assert( pNode->pCutBest[0] || pNode->pCutBest[1] ); // at least one assigned
//    return pNode->nRefEst[0] + pNode->nRefEst[1];
    return pNode->nRefEst[2];
}


/**Function*************************************************************

  Synopsis    [Increments the actual reference counter of a phase.]

  Description [Returns the old reference counter.]
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int Map_NodeIncRefPhaseAct( Map_Node_t * pNode, int fPhase )
{
    assert( !Map_IsComplement(pNode) );
    if ( pNode->pCutBest[0] && pNode->pCutBest[1] ) // both assigned
        return pNode->nRefAct[fPhase]++;
    assert( pNode->pCutBest[0] || pNode->pCutBest[1] ); // at least one assigned
    return pNode->nRefAct[2]++;
}

/**Function*************************************************************

  Synopsis    [Decrements the actual reference counter of a phase.]

  Description [Returns the new reference counter.]
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int Map_NodeDecRefPhaseAct( Map_Node_t * pNode, int fPhase )
{
    assert( !Map_IsComplement(pNode) );
    if ( pNode->pCutBest[0] && pNode->pCutBest[1] ) // both assigned
        return --pNode->nRefAct[fPhase];
    assert( pNode->pCutBest[0] || pNode->pCutBest[1] ); // at least one assigned
    return --pNode->nRefAct[2];
}


/**Function*************************************************************

  Synopsis    [Sets the estimated reference counter for the PIs.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Map_MappingEstimateRefsInit( Map_Man_t * p )
{
    Map_Node_t * pNode;
    int i;
    for ( i = 0; i < p->vMapObjs->nSize; i++ )
    {
        pNode = p->vMapObjs->pArray[i];
//        pNode->nRefEst[0] = pNode->nRefEst[1] = ((float)pNode->nRefs)*(float)2.0;
        pNode->nRefEst[0] = pNode->nRefEst[1] = pNode->nRefEst[2] = ((float)pNode->nRefs);
    }
}

/**Function*************************************************************

  Synopsis    [Sets the estimated reference counter.]

  Description [When this procedure is called for the first time,
  the reference counter is estimated from the AIG. Otherwise, it is
  a linear combination of reference counters in the last two iterations.]
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Map_MappingEstimateRefs( Map_Man_t * p )
{
    Map_Node_t * pNode;
    int i;
    for ( i = 0; i < p->vMapObjs->nSize; i++ )
    {
        pNode = p->vMapObjs->pArray[i];
//        pNode->nRefEst[0] = (float)((2.0 * pNode->nRefEst[0] + 1.0 * pNode->nRefAct[0]) / 3.0);
//        pNode->nRefEst[1] = (float)((2.0 * pNode->nRefEst[1] + 1.0 * pNode->nRefAct[1]) / 3.0);
//        pNode->nRefEst[2] = (float)((2.0 * pNode->nRefEst[2] + 1.0 * pNode->nRefAct[2]) / 3.0);
        pNode->nRefEst[0] = (float)((3.0 * pNode->nRefEst[0] + 1.0 * pNode->nRefAct[0]) / 4.0);
        pNode->nRefEst[1] = (float)((3.0 * pNode->nRefEst[1] + 1.0 * pNode->nRefAct[1]) / 4.0);
        pNode->nRefEst[2] = (float)((3.0 * pNode->nRefEst[2] + 1.0 * pNode->nRefAct[2]) / 4.0);
    }
}

/**function*************************************************************

  synopsis    [Computes the area flow of the cut.]

  description [Computes the area flow of the cut if it is implemented using 
  the best supergate with the best phase.]
               
  sideeffects []

  seealso     []

***********************************************************************/
float Map_CutGetAreaFlow( Map_Cut_t * pCut, int fPhase )
{
    Map_Match_t * pM = pCut->M + fPhase;
    Map_Super_t * pSuper = pM->pSuperBest;
    unsigned uPhaseTot = pM->uPhaseBest;
    Map_Cut_t * pCutFanin;
    float aFlowRes, aFlowFanin, nRefs;
    int i, fPinPhasePos;

    // start the resulting area flow
    aFlowRes = pSuper->Area;
    // iterate through the leaves
    for ( i = 0; i < pCut->nLeaves; i++ )
    {
        // get the phase of this fanin
        fPinPhasePos = ((uPhaseTot & (1 << i)) == 0);
        // get the cut implementing this phase of the fanin
        pCutFanin = pCut->ppLeaves[i]->pCutBest[fPinPhasePos];
        // if the cut is not available, we have to use the opposite phase
        if ( pCutFanin == NULL )
        {
            fPinPhasePos = !fPinPhasePos;
            pCutFanin = pCut->ppLeaves[i]->pCutBest[fPinPhasePos];
        }
        aFlowFanin = pCutFanin->M[fPinPhasePos].AreaFlow; // ignores the area of the interter
        // get the fanout count of the cut in the given phase
        nRefs = Map_NodeReadRefPhaseEst( pCut->ppLeaves[i], fPinPhasePos );
        // if the node does no fanout, assume fanout count equal to 1
        if ( nRefs == (float)0.0 )
            nRefs = (float)1.0;
        // add the area flow due to the fanin
        aFlowRes += aFlowFanin / nRefs;
    }
    pM->AreaFlow = aFlowRes;
    return aFlowRes;
}

/**function*************************************************************

  synopsis    [Computes the power flow of the cut.]

  description [Computes the power flow of the cut if it is implemented using
  the best supergate with the best phase.]

  sideeffects []

  seealso     []

***********************************************************************/
float Map_CutGetPowerFlow(  Map_Node_t * pNode, Map_Cut_t * pCut, int fPhase )
{
    Map_Match_t * pM = pCut->M + fPhase;
    Map_Super_t * pSuper = pM->pSuperBest;
    Map_Node_t * pNodeChild;
    unsigned uPhaseTot = pM->uPhaseBest;
    Map_Cut_t * pCutFanin;
    float aFlowRes, aFlowFanin, nRefs, aSwitchActivityInt = 1;
    int i, fPinPhasePos;

    // start the power of this cut
    // net switching power of the root gate
    aFlowRes = pSuper->PowerSwi * pNode->Switching;

    // Add cell internal power (dependent on internal power of the root gate and the switching activity seen at the input pins)
    // Use Formula: aSwitchingInternal = powerInt ( 1 - pNode->Switching ) * ( 1 - ( 1 - pNodeChild1->Switching ) * ( pNodeChild2->Switching ) )

    for ( i = 0; i < pCut->nLeaves; i++ )
    {
        pNodeChild  = pCut->ppLeaves[i];
        // For cell internal switching: Calculate the probability of no child node switching
        aSwitchActivityInt *= 1 - pNodeChild->Switching;
    }

    aFlowRes += pSuper->PowerInt * ( 1 - pNode->Switching ) * ( 1 - aSwitchActivityInt );
    // iterate through the leaves
    for ( i = 0; i < pCut->nLeaves; i++ )
    {
        // get the phase of this fanin
        fPinPhasePos = ((uPhaseTot & (1 << i)) == 0);
        // get the cut implementing this phase of the fanin
        pCutFanin = pCut->ppLeaves[i]->pCutBest[fPinPhasePos];
        // if the cut is not available, we have to use the opposite phase
        if ( pCutFanin == NULL )
        {
            fPinPhasePos = !fPinPhasePos;
            pCutFanin = pCut->ppLeaves[i]->pCutBest[fPinPhasePos];
        }
        aFlowFanin = pCutFanin->M[fPinPhasePos].PowerF; // ignores the power of the inverter
        // get the fanout count of the cut in the given phase
        nRefs = Map_NodeReadRefPhaseEst( pCut->ppLeaves[i], fPinPhasePos );
        // if the node does no fanout, assume fanout count equal to 1
        if ( nRefs == (float)0.0 )
            nRefs = (float)1.0;
        // add the power due to the fanin
        aFlowRes += aFlowFanin / nRefs;
    }
    pM->PowerF = aFlowRes;
    return aFlowRes;
}

/**function*************************************************************

  synopsis    [References or dereferences the cut.]

  description [This reference part is similar to Cudd_NodeReclaim(). 
  The dereference part is similar to Cudd_RecursiveDeref().]
               
  sideeffects []

  seealso     []

***********************************************************************/
float Map_CutRefDeref( Map_Cut_t * pCut, int fPhase, int fReference, int fUpdateProf )
{
    Map_Node_t * pNodeChild;
    Map_Cut_t * pCutChild;
    float aArea;
    int i, fPhaseChild;
//    int nRefs;

    // consider the elementary variable
    if ( pCut->nLeaves == 1 )
        return 0;
    // start the area of this cut
    aArea = Map_CutGetRootArea( pCut, fPhase );
    if ( fUpdateProf )
    {
        if ( fReference )
            Mio_GateIncProfile2( pCut->M[fPhase].pSuperBest->pRoot );
        else
            Mio_GateDecProfile2( pCut->M[fPhase].pSuperBest->pRoot );
    }
    // go through the children
    for ( i = 0; i < pCut->nLeaves; i++ )
    {
        pNodeChild  = pCut->ppLeaves[i];
        fPhaseChild = Map_CutGetLeafPhase( pCut, fPhase, i );
        // get the reference counter of the child
/*
        // this code does not take inverters into account
        // the quality of area recovery seems to always be a little worse
        if ( fReference )
            nRefs = Map_NodeIncRefPhaseAct( pNodeChild, fPhaseChild );
        else
            nRefs = Map_NodeDecRefPhaseAct( pNodeChild, fPhaseChild );
        assert( nRefs >= 0 );
        // skip if the child was already reference before
        if ( nRefs > 0 )
            continue;
*/

        if ( fReference )
        {
            if ( pNodeChild->pCutBest[0] && pNodeChild->pCutBest[1] ) // both phases are present
            {
                // if this phase of the node is referenced, there is no recursive call
                pNodeChild->nRefAct[2]++;
                if ( pNodeChild->nRefAct[fPhaseChild]++ > 0 )
                    continue;
            }
            else // only one phase is present
            {
                // inverter should be added if the phase
                // (a) has no reference and (b) is implemented using other phase
                if ( pNodeChild->nRefAct[fPhaseChild]++ == 0 && pNodeChild->pCutBest[fPhaseChild] == NULL )
                    aArea += pNodeChild->p->pSuperLib->AreaInv;
                // if the node is referenced, there is no recursive call
                if ( pNodeChild->nRefAct[2]++ > 0 )
                    continue;
            }
        }
        else
        {
            if ( pNodeChild->pCutBest[0] && pNodeChild->pCutBest[1] ) // both phases are present
            {
                // if this phase of the node is referenced, there is no recursive call
                --pNodeChild->nRefAct[2];
                if ( --pNodeChild->nRefAct[fPhaseChild] > 0 )
                    continue;
            }
            else // only one phase is present
            {
                // inverter should be added if the phase
                // (a) has no reference and (b) is implemented using other phase
                if ( --pNodeChild->nRefAct[fPhaseChild] == 0 && pNodeChild->pCutBest[fPhaseChild] == NULL )
                    aArea += pNodeChild->p->pSuperLib->AreaInv;
                // if the node is referenced, there is no recursive call
                if ( --pNodeChild->nRefAct[2] > 0 )
                    continue;
            }
            assert( pNodeChild->nRefAct[fPhaseChild] >= 0 );
        }

        // get the child cut
        pCutChild = pNodeChild->pCutBest[fPhaseChild];
        // if the child does not have this phase mapped, take the opposite phase
        if ( pCutChild == NULL )
        {
            fPhaseChild = !fPhaseChild;
            pCutChild   = pNodeChild->pCutBest[fPhaseChild];
        }
        // reference and compute area recursively
        aArea += Map_CutRefDeref( pCutChild, fPhaseChild, fReference, fUpdateProf );
    }
    return aArea;
}

/**function*************************************************************

  synopsis    [References or dereferences the cut.]

  description [This reference part is similar to Cudd_NodeReclaim().
  The dereference part is similar to Cudd_RecursiveDeref().]

  sideeffects []

  seealso     []

***********************************************************************/
float Map_PowerCutRefDeref( Map_Node_t * pNode, Map_Cut_t * pCut, int fPhase, int fReference, int fUpdateProf )
{
    Map_Node_t * pNodeChild;
    Map_Cut_t * pCutChild;
    float aPower, aSwitchActivity, aSwitchActivityInt = 1;
    int i, fPhaseChild;
//    int nRefs;
    // start switching activity for the node
    aSwitchActivity = pNode->Switching;
    // consider the elementary variable
    if ( pCut->nLeaves == 1 )
        return aSwitchActivity;
    // start the power of this cut
    // net switching power of the root gate
    aPower = Map_CutGetRootPowerSwi( pCut, fPhase ) * aSwitchActivity;

    // Add cell internal power (dependent on internal power of the root gate and the switching activity seen at the input pins)
    // Use Formula: aSwitchingInternal = powerInt ( 1 - pNode->Switching ) * ( 1 - ( 1 - pNodeChild1->Switching ) * ( pNodeChild2->Switching ) )

    for ( i = 0; i < pCut->nLeaves; i++ )
    {
        pNodeChild  = pCut->ppLeaves[i];
        // For cell internal switching: Calculate the probability of no child node switching
        aSwitchActivityInt *= 1 - pNodeChild->Switching;
    }
    // Use Formula: Switching Internal = 1 - P(output switching) - P(no child node switching)
    aSwitchActivityInt = ( 1 - aSwitchActivity - aSwitchActivityInt );
    aPower += Map_CutGetRootPowerInt( pCut, fPhase ) * aSwitchActivityInt;

    if ( fUpdateProf )
    {
        if ( fReference )
            Mio_GateIncProfile2( pCut->M[fPhase].pSuperBest->pRoot );
        else
            Mio_GateDecProfile2( pCut->M[fPhase].pSuperBest->pRoot );
    }
    // go through the children
    for ( i = 0; i < pCut->nLeaves; i++ )
    {
        pNodeChild  = pCut->ppLeaves[i];
        fPhaseChild = Map_CutGetLeafPhase( pCut, fPhase, i );

        //Add Power dependent on Phase
        if ( fReference )
        {
            if ( pNodeChild->pCutBest[0] && pNodeChild->pCutBest[1] ) // both phases are present
            {
                // if this phase of the node is referenced, there is no recursive call
                pNodeChild->nRefAct[2]++;
                if ( pNodeChild->nRefAct[fPhaseChild]++ > 0 )
                    continue;
            }
            else // only one phase is present
            {
                // inverter should be added if the phase
                // (a) has no reference and (b) is implemented using other phase
                // Inverter only dissipates net switching power
                if ( pNodeChild->nRefAct[fPhaseChild]++ == 0 && pNodeChild->pCutBest[fPhaseChild] == NULL )
                    aPower += pNodeChild->p->pSuperLib->PowerInv * pNodeChild->Switching;
                // if the node is referenced, there is no recursive call
                if ( pNodeChild->nRefAct[2]++ > 0 )
                    continue;
            }
        }
        else
        {
            if ( pNodeChild->pCutBest[0] && pNodeChild->pCutBest[1] ) // both phases are present
            {
                // if this phase of the node is referenced, there is no recursive call
                --pNodeChild->nRefAct[2];
                if ( --pNodeChild->nRefAct[fPhaseChild] > 0 )
                    continue;
            }
            else // only one phase is present
            {
                // inverter should be added if the phase
                // (a) has no reference and (b) is implemented using other phase
                // todo: when splitting up the power the inverter needs to include cellint and netswi power
                if ( --pNodeChild->nRefAct[fPhaseChild] == 0 && pNodeChild->pCutBest[fPhaseChild] == NULL )
                    aPower += pNodeChild->p->pSuperLib->PowerInv * pNodeChild->Switching;
                // if the node is referenced, there is no recursive call
                if ( --pNodeChild->nRefAct[2] > 0 )
                    continue;
            }
            assert( pNodeChild->nRefAct[fPhaseChild] >= 0 );
        }

        // get the child cut
        pCutChild = pNodeChild->pCutBest[fPhaseChild];
        // if the child does not have this phase mapped, take the opposite phase
        if ( pCutChild == NULL )
        {
            fPhaseChild = !fPhaseChild;
            pCutChild   = pNodeChild->pCutBest[fPhaseChild];
        }

        // reference and compute power recursively
        aPower += Map_PowerCutRefDeref( pNodeChild, pCutChild, fPhaseChild, fReference, fUpdateProf );
    }
    return aPower;
}

/**function*************************************************************

  synopsis    [Computes the exact area associated with the cut.]

  description [Assumes that the cut is referenced.]
               
  sideeffects []

  seealso     []

***********************************************************************/
float Map_CutGetAreaRefed( Map_Cut_t * pCut, int fPhase )
{
    float aResult, aResult2;
    aResult2 = Map_CutRefDeref( pCut, fPhase, 0, 0 ); // dereference
    aResult  = Map_CutRefDeref( pCut, fPhase, 1, 0 ); // reference
//    assert( aResult == aResult2 );
    return aResult;
}

/**function*************************************************************

  synopsis    [Computes the exact area associated with the cut.]

  description []
               
  sideeffects []

  seealso     []

***********************************************************************/
float Map_CutGetAreaDerefed( Map_Cut_t * pCut, int fPhase )
{
    float aResult, aResult2;
    aResult2 = Map_CutRefDeref( pCut, fPhase, 1, 0 ); // reference
    aResult  = Map_CutRefDeref( pCut, fPhase, 0, 0 ); // dereference
//    assert( aResult == aResult2 );
    return aResult;
}

/**function*************************************************************

  synopsis    [Computes the exact power associated with the cut.]

  description []

  sideeffects []

  seealso     []

***********************************************************************/
float Map_CutGetPowerDerefed( Map_Node_t * pNode, Map_Cut_t * pCut, int fPhase )
{
    float aResult, aResult2;
    aResult2 = Map_PowerCutRefDeref( pNode, pCut, fPhase, 1, 0 ); // reference
    aResult  = Map_PowerCutRefDeref( pNode, pCut, fPhase, 0, 0 ); // dereference
//    assert( aResult == aResult2 );
    return aResult;
}

/**function*************************************************************

  synopsis    [References the cut.]

  description []
               
  sideeffects []

  seealso     []

***********************************************************************/
float Map_CutRef( Map_Cut_t * pCut, int fPhase, int fProfile )
{
    return Map_CutRefDeref( pCut, fPhase, 1, fProfile ); // reference
}

/**function*************************************************************

  synopsis    [Dereferences the cut.]

  description []
               
  sideeffects []

  seealso     []

***********************************************************************/
float Map_CutDeref( Map_Cut_t * pCut, int fPhase, int fProfile )
{
    return Map_CutRefDeref( pCut, fPhase, 0, fProfile ); // dereference
}

/**function*************************************************************

  synopsis    [References the cut.]

  description []

  sideeffects []

  seealso     []

***********************************************************************/
float Map_PowerCutRef( Map_Node_t * pNode, Map_Cut_t * pCut, int fPhase, int fProfile )
{
    return Map_PowerCutRefDeref( pNode, pCut, fPhase, 1, fProfile ); // reference
}

/**function*************************************************************

  synopsis    [Dereferences the cut.]

  description []

  sideeffects []

  seealso     []

***********************************************************************/
float Map_PowerCutDeref( Map_Node_t * pNode, Map_Cut_t * pCut, int fPhase, int fProfile )
{
    return Map_PowerCutRefDeref( pNode, pCut, fPhase, 0, fProfile ); // dereference
}

/**Function*************************************************************

  Synopsis    [Computes actual reference counters.]

  Description [Collects the nodes used in the mapping in array pMan->vMapping.
  Nodes are collected in reverse topological order to facilitate the 
  computation of required times.]
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Map_MappingSetRefs_rec( Map_Man_t * pMan, Map_Node_t * pNode )
{
    Map_Cut_t * pCut;
    Map_Node_t * pNodeR;
    unsigned uPhase;
    int i, fPhase, fInvPin;
    // get the regular node and its phase
    pNodeR = Map_Regular(pNode);
    fPhase = !Map_IsComplement(pNode);
    pNodeR->nRefAct[2]++;
    // quit if the node was already visited in this phase
    if ( pNodeR->nRefAct[fPhase]++ )
        return;
    // quit if this is a PI node
    if ( Map_NodeIsVar(pNodeR) )
        return;
    // propagate through buffer
    if ( Map_NodeIsBuf(pNodeR) )
    {
        Map_MappingSetRefs_rec( pMan, Map_NotCond(pNodeR->p1, Map_IsComplement(pNode)) );
        return;
    }
    assert( Map_NodeIsAnd(pNode) );
    // get the cut implementing this or opposite polarity
    pCut = pNodeR->pCutBest[fPhase];
    if ( pCut == NULL )
    {
        fPhase = !fPhase;
        pCut   = pNodeR->pCutBest[fPhase];
    }
    if ( pMan->fUseProfile )
        Mio_GateIncProfile2( pCut->M[fPhase].pSuperBest->pRoot );
    // visit the transitive fanin
    uPhase = pCut->M[fPhase].uPhaseBest;
    for ( i = 0; i < pCut->nLeaves; i++ )
    {
        fInvPin = ((uPhase & (1 << i)) > 0);
        Map_MappingSetRefs_rec( pMan, Map_NotCond(pCut->ppLeaves[i], fInvPin) );
    }
}
void Map_MappingSetRefs( Map_Man_t * pMan )
{
    Map_Node_t * pNode;
    int i;
    if ( pMan->fUseProfile )
        Mio_LibraryCleanProfile2( pMan->pSuperLib->pGenlib );
    // clean all references
    for ( i = 0; i < pMan->vMapObjs->nSize; i++ )
    {
        pNode = pMan->vMapObjs->pArray[i];
        pNode->nRefAct[0] = 0;
        pNode->nRefAct[1] = 0;
        pNode->nRefAct[2] = 0;
    }
    // visit nodes reachable from POs in the DFS order through the best cuts
    for ( i = 0; i < pMan->nOutputs; i++ )
    {
        pNode = pMan->pOutputs[i];
        if ( !Map_NodeIsConst(pNode) )
            Map_MappingSetRefs_rec( pMan, pNode );
    }
}


/**Function*************************************************************

  Synopsis    [Computes the array of mapping.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
float Map_MappingGetArea( Map_Man_t * pMan )
{
    Map_Node_t * pNode;
    float Area = 0.0;
    int i;
    if ( pMan->fUseProfile )
        Mio_LibraryCleanProfile2( pMan->pSuperLib->pGenlib );
    for ( i = 0; i < pMan->vMapObjs->nSize; i++ )
    {
        pNode = pMan->vMapObjs->pArray[i];
        if ( pNode->nRefAct[2] == 0 )
            continue;
        if ( Map_NodeIsBuf(pNode) )
            continue;
        // at least one phase has the best cut assigned
        assert( pNode->pCutBest[0] != NULL || pNode->pCutBest[1] != NULL );
        // at least one phase is used in the mapping
        assert( pNode->nRefAct[0] > 0 || pNode->nRefAct[1] > 0 );
        // compute the array due to the supergate
        if ( Map_NodeIsAnd(pNode) )
        {
            // count area of the negative phase
            if ( pNode->pCutBest[0] && (pNode->nRefAct[0] > 0 || pNode->pCutBest[1] == NULL) )
            {
                Area += pNode->pCutBest[0]->M[0].pSuperBest->Area;
                if ( pMan->fUseProfile )
                    Mio_GateIncProfile2( pNode->pCutBest[0]->M[0].pSuperBest->pRoot );
            }
            // count area of the positive phase
            if ( pNode->pCutBest[1] && (pNode->nRefAct[1] > 0 || pNode->pCutBest[0] == NULL) )
            {
                Area += pNode->pCutBest[1]->M[1].pSuperBest->Area;
                if ( pMan->fUseProfile )
                    Mio_GateIncProfile2( pNode->pCutBest[1]->M[1].pSuperBest->pRoot );
            }
        }
        // count area of the interver if we need to implement one phase with another phase
        if ( (pNode->pCutBest[0] == NULL && pNode->nRefAct[0] > 0) || 
             (pNode->pCutBest[1] == NULL && pNode->nRefAct[1] > 0) )
            Area += pMan->pSuperLib->AreaInv;
    }
    // add buffers for each CO driven by a CI
    for ( i = 0; i < pMan->nOutputs; i++ )
        if ( Map_NodeIsVar(pMan->pOutputs[i]) && !Map_IsComplement(pMan->pOutputs[i]) )
            Area += pMan->pSuperLib->AreaBuf;
    return Area;
}

/**Function*************************************************************

  Synopsis    [Computes the power of mapping.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
float Map_MappingGetPower( Map_Man_t * pMan )
{
    Map_Node_t * pNode, * pNodeChild;
    float Power = 0.0, aSwitchActivityInt = 1.0;
    int i, j;
    if ( pMan->fUseProfile )
        Mio_LibraryCleanProfile2( pMan->pSuperLib->pGenlib );
    for ( i = 0; i < pMan->vMapObjs->nSize; i++ )
    {
        pNode = pMan->vMapObjs->pArray[i];
        if ( pNode->nRefAct[2] == 0 )
            continue;
        if ( Map_NodeIsBuf(pNode) )
            continue;
        // at least one phase has the best cut assigned
        assert( pNode->pCutBest[0] != NULL || pNode->pCutBest[1] != NULL );
        // at least one phase is used in the mapping
        assert( pNode->nRefAct[0] > 0 || pNode->nRefAct[1] > 0 );
        // compute the array due to the supergate
        aSwitchActivityInt = 1.0;
        if ( Map_NodeIsAnd(pNode) )
        {
            // count power of the negative phase
            if ( pNode->pCutBest[0] && (pNode->nRefAct[0] > 0 || pNode->pCutBest[1] == NULL) )
            {
                Power += pNode->pCutBest[0]->M[0].pSuperBest->PowerSwi * pNode->Switching;
                for ( j = 0; j < pNode->pCutBest[0]->nLeaves; j++ )
                {
                    pNodeChild  = pNode->pCutBest[0]->ppLeaves[j];
                    // For cell internal switching: Calculate the probability of no child node switching
                    aSwitchActivityInt *= 1 - pNodeChild->Switching;
                }
                // Use Formula: Switching Internal = 1 - P(output switching) - P(no child node switching)
                aSwitchActivityInt = 1 - pNode->Switching - aSwitchActivityInt;
                Power += pNode->pCutBest[0]->M[0].pSuperBest->PowerInt * aSwitchActivityInt;

                if ( pMan->fUseProfile )
                    Mio_GateIncProfile2( pNode->pCutBest[0]->M[0].pSuperBest->pRoot );
            }
            // count power of the positive phase
            if ( pNode->pCutBest[1] && (pNode->nRefAct[1] > 0 || pNode->pCutBest[0] == NULL) )
            {
                Power += pNode->pCutBest[1]->M[1].pSuperBest->PowerSwi * pNode->Switching;
                for ( j = 0; j < pNode->pCutBest[1]->nLeaves; j++ )
                {
                    pNodeChild  = pNode->pCutBest[1]->ppLeaves[j];
                    // For cell internal switching: Calculate the probability of no child node switching
                    aSwitchActivityInt *= 1 - pNodeChild->Switching;
                }
                // Use Formula: Switching Internal = 1 - P(output switching) - P(no child node switching)
                aSwitchActivityInt = 1 - pNode->Switching - aSwitchActivityInt;
                Power += pNode->pCutBest[1]->M[1].pSuperBest->PowerInt * aSwitchActivityInt;

                if ( pMan->fUseProfile )
                    Mio_GateIncProfile2( pNode->pCutBest[1]->M[1].pSuperBest->pRoot );
            }
        }
        // count power of the inverter if we need to implement one phase with another phase
        if ( (pNode->pCutBest[0] == NULL && pNode->nRefAct[0] > 0) ||
             (pNode->pCutBest[1] == NULL && pNode->nRefAct[1] > 0) )
            Power += pMan->pSuperLib->PowerInv * pNode->Switching;
    }
    // add buffers for each CO driven by a CI
    for ( i = 0; i < pMan->nOutputs; i++ )
        if ( Map_NodeIsVar(pMan->pOutputs[i]) && !Map_IsComplement(pMan->pOutputs[i]) )
            Power += pMan->pSuperLib->PowerBuf * pNode->Switching;
    return Power;
}

/**Function*************************************************************

  Synopsis    [Computes the power of mapping.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
    /*ToDO(hienbenj): get the net switching power of a gate and if it is high look at the fanins of its fanout gate*/
float Map_MappingGetSwitchingDiff( Map_Man_t * pMan )
{
    Mio_Gate_t * Gate;
    Map_Node_t * pNode, * pNodeChild;
    float Power = 0.0, aSwitchActivityInt = 1.0, PowerInt = 0.0, PowerSwi = 0.0;
    int i, j;
    if ( pMan->fUseProfile )
        Mio_LibraryCleanProfile2( pMan->pSuperLib->pGenlib );
    for ( i = 0; i < pMan->vMapObjs->nSize; i++ )
    {
        pNode = pMan->vMapObjs->pArray[i];
        if ( pNode->nRefAct[2] == 0 )
            continue;
        if ( Map_NodeIsBuf(pNode) )
            continue;
        // at least one phase has the best cut assigned
        assert( pNode->pCutBest[0] != NULL || pNode->pCutBest[1] != NULL );
        // at least one phase is used in the mapping
        assert( pNode->nRefAct[0] > 0 || pNode->nRefAct[1] > 0 );
        // compute the array due to the supergate
        aSwitchActivityInt = 1.0;
        if ( Map_NodeIsAnd(pNode) )
        {
            // count power of the negative phase
            if ( pNode->pCutBest[0] && (pNode->nRefAct[0] > 0 || pNode->pCutBest[1] == NULL) )
            {
                PowerSwi = pNode->pCutBest[0]->M[0].pSuperBest->PowerSwi;
                Power += PowerSwi * pNode->Switching;
                for ( j = 0; j < pNode->pCutBest[0]->nLeaves; j++ )
                {
                    pNodeChild  = pNode->pCutBest[0]->ppLeaves[j];
                    // For cell internal switching: Calculate the probability of no child node switching
                    aSwitchActivityInt = pNodeChild->Switching;
                }
                // Use Formula: Switching Internal = 1 - P(output switching) - P(no child node switching)
                aSwitchActivityInt = 1 - pNode->Switching - aSwitchActivityInt;
                PowerInt = pNode->pCutBest[0]->M[0].pSuperBest->PowerInt * aSwitchActivityInt;
                Power += PowerInt;

                if ( pMan->fUseProfile )
                    Mio_GateIncProfile2( pNode->pCutBest[0]->M[0].pSuperBest->pRoot );

                if(PowerSwi > 0.5)
                {
                    Gate = pNode->pCutBest[0]->M[0].pSuperBest->pRoot;
                    volatile int i = 0;
                }
            }
            // count power of the positive phase
            if ( pNode->pCutBest[1] && (pNode->nRefAct[1] > 0 || pNode->pCutBest[0] == NULL) )
            {
                PowerSwi = pNode->pCutBest[1]->M[1].pSuperBest->PowerSwi;
                Power += PowerSwi * pNode->Switching;
                for ( j = 0; j < pNode->pCutBest[1]->nLeaves; j++ )
                {
                    pNodeChild  = pNode->pCutBest[1]->ppLeaves[j];
                    // For cell internal switching: Calculate the probability of no child node switching
                    aSwitchActivityInt *= 1 - pNodeChild->Switching;
                }
                // Use Formula: Switching Internal = 1 - P(output switching) - P(no child node switching)
                aSwitchActivityInt = 1 - pNode->Switching - aSwitchActivityInt;
                PowerInt = pNode->pCutBest[1]->M[1].pSuperBest->PowerInt * aSwitchActivityInt;
                Power += PowerInt;

                if ( pMan->fUseProfile )
                    Mio_GateIncProfile2( pNode->pCutBest[1]->M[1].pSuperBest->pRoot );

                if(PowerSwi > 0.5)
                {
                    Gate = pNode->pCutBest[1]->M[1].pSuperBest->pRoot;
                    volatile int i = 0;
                }
            }
        }
    }
    return Power;
}

////////////////////////////////////////////////////////////////////////
///                       END OF FILE                                ///
////////////////////////////////////////////////////////////////////////


ABC_NAMESPACE_IMPL_END

