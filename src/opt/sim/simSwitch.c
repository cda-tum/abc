/**CFile****************************************************************

  FileName    [simSwitch.c]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [Network and node package.]

  Synopsis    [Computes switching activity of nodes in the ABC network.]

  Author      [Alan Mishchenko]
  
  Affiliation [UC Berkeley]

  Date        [Ver. 1.0. Started - June 20, 2005.]

  Revision    [$Id: simSwitch.c,v 1.00 2005/06/20 00:00:00 alanmi Exp $]

***********************************************************************/

#include <math.h>
#include "base/abc/abc.h"
#include "sim.h"

ABC_NAMESPACE_IMPL_START


////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

static void Sim_NodeSimulate( Abc_Obj_t * pNode, Vec_Ptr_t * vSimInfo, int nSimWords );
static float Sim_ComputeSwitching( unsigned * pSimInfo, int nSimWords );
static float Sim_ComputeSwitchingT( Abc_Obj_t * pNode, Vec_Ptr_t * vSimInfo, int nSimWords );

////////////////////////////////////////////////////////////////////////
///                     FUNCTION DEFINITIONS                         ///
////////////////////////////////////////////////////////////////////////

/**Function*************************************************************

  Synopsis    [Computes switching activity using simulation.]

  Description [Computes switching activity, which is understood as the
  probability of switching under random simulation. Assigns the
  random simulation information at the CI and propagates it through
  the internal nodes of the AIG.]
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
Vec_Int_t * Sim_NtkComputeSwitching( Abc_Ntk_t * pNtk, int nPatterns )
{
    Vec_Int_t * vSwitching;
    float * pSwitching;
    Vec_Int_t * vIntSwitching;
    float * pIntSwitching;
    Vec_Ptr_t * vNodes;
    Vec_Ptr_t * vSimInfo;
    Abc_Obj_t * pNode;
    unsigned * pSimInfo;
    int nSimWords, i;
    float output = 0;

    // allocate space for simulation info of all nodes
    nSimWords = SIM_NUM_WORDS(nPatterns);
    vSimInfo = Sim_UtilInfoAlloc( Abc_NtkObjNumMax(pNtk), nSimWords, 0 );
    // assign the random simulation to the CIs
    vSwitching = Vec_IntStart( Abc_NtkObjNumMax(pNtk) );
    pSwitching = (float *)vSwitching->pArray;
    vIntSwitching = Vec_IntStart( Abc_NtkObjNumMax(pNtk) );
    pIntSwitching = (float *)vIntSwitching->pArray;
    Abc_NtkForEachCi( pNtk, pNode, i )
    {
        //printf("I\n");
        pSimInfo = (unsigned *)Vec_PtrEntry(vSimInfo, pNode->Id);
        /*if(pNode->Id == 1)
        {
            Sim_UtilSetBiased( pSimInfo, nSimWords );
        }
        else
        {*/
            Sim_UtilSetRandom( pSimInfo, nSimWords );
        /*}*/

        pSwitching[pNode->Id] = Sim_ComputeSwitching( pSimInfo, nSimWords );
        output = pSwitching[pNode->Id];
    }
    // simulate the internal nodes
    vNodes  = Abc_AigDfs( pNtk, 1, 0 );
    Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pNode, i )
    {
        //printf("N\n");
        pSimInfo = (unsigned *)Vec_PtrEntry(vSimInfo, pNode->Id);
        Sim_UtilSimulateNodeOne( pNode, vSimInfo, nSimWords, 0 );
        pIntSwitching[pNode->Id] = Sim_ComputeSwitchingT( pNode, vSimInfo, nSimWords );
        output = (pSwitching[Abc_ObjFaninId1(pNode)] - pSwitching[Abc_ObjFaninId0(pNode)]) / pSwitching[Abc_ObjFaninId1(pNode)];
        printf("Percentage Diff %f\n", output*100);
        output = pIntSwitching[pNode->Id];
        printf("Transitions %f\n", output);
        pSwitching[pNode->Id] = Sim_ComputeSwitching( pSimInfo, nSimWords );
        output = pSwitching[pNode->Id];
        printf("Switching %f\n", output);
    }
    Vec_PtrFree( vNodes );
    Sim_UtilInfoFree( vSimInfo );
    return vSwitching;
}

/**Function*************************************************************

  Synopsis    [Computes switching activity of one node.]

  Description [Uses the formula: Switching = 2 * nOnes * nZeros / (nTotal ^ 2) ]
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
float Sim_ComputeSwitching( unsigned * pSimInfo, int nSimWords )
{
    int nOnes, nTotal;
    nTotal = 32 * nSimWords;
    nOnes = Sim_UtilCountOnes( pSimInfo, nSimWords );
    return (float)2.0 * nOnes / nTotal * (nTotal - nOnes) / nTotal;
}

/**Function*************************************************************

  Synopsis    [Computes switching activity of one node.]

  Description [Uses the formula: Switching = 2 * nOnes * nZeros / (nTotal ^ 2) ]

  SideEffects []

  SeeAlso     []

***********************************************************************/
float Sim_ComputeSwitchingT( Abc_Obj_t * pNode, Vec_Ptr_t * vSimInfo, int nSimWords )
{
    int nOnes, nTotal, nTransitions;
    float Swi, Swi_ones;
    unsigned * pSimmNode, * pSimmNode1, * pSimmNode2;

    pSimmNode  = (unsigned *)Vec_PtrEntry(vSimInfo, pNode->Id);
    pSimmNode1 = (unsigned *)Vec_PtrEntry(vSimInfo, Abc_ObjFaninId0(pNode));
    pSimmNode2 = (unsigned *)Vec_PtrEntry(vSimInfo, Abc_ObjFaninId1(pNode));

    nTotal = 32 * nSimWords;
    nOnes = Sim_UtilCountOnes( pSimmNode, nSimWords );

    nTransitions = Sim_UtilCountTransitions( pSimmNode, nSimWords );
    Swi = (float) nTransitions/nTotal;
    //printf("Switching%f\n", Swi);
    nTransitions = Sim_UtilCountCellIntTransitions( pSimmNode, pSimmNode1, pSimmNode2, nSimWords );
    Swi = (float) nTransitions/nTotal;
    //printf("CellInt %f\n", Swi);
    /*nTotal = 32 * nSimWords;
    Swi_ones = (float)2.0 * nOnes / nTotal * (nTotal - nOnes) / nTotal;
    printf("%f\n", Swi_ones);
    Swi = Swi - Swi_ones;
    printf("%f\n", Swi);*/
    //round(x * 8.0) / 8.0;

    return Swi;
}

////////////////////////////////////////////////////////////////////////
///                       END OF FILE                                ///
////////////////////////////////////////////////////////////////////////


ABC_NAMESPACE_IMPL_END

