/**CFile****************************************************************

  FileName    [abcIf.c]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [Network and node package.]

  Synopsis    [Interface with the FPGA mapping package.]

  Author      [Alan Mishchenko]
  
  Affiliation [UC Berkeley]

  Date        [Ver. 1.0. Started - November 21, 2006.]

  Revision    [$Id: abcIf.c,v 1.00 2006/11/21 00:00:00 alanmi Exp $]

***********************************************************************/

#include "base/abc/abc.h"
#include "base/main/main.h"
#include "map/if/if.h"
#include "bool/kit/kit.h"
#include "aig/aig/aig.h"
#include "map/mio/mio.h"

ABC_NAMESPACE_IMPL_START


////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

extern If_Man_t *  Abc_NtkToIf( Abc_Ntk_t * pNtk, If_Par_t * pPars );
static Abc_Ntk_t * Abc_NtkFromIf( If_Man_t * pIfMan, Abc_Ntk_t * pNtk );
extern Abc_Obj_t * Abc_NodeFromIf_rec( Abc_Ntk_t * pNtkNew, If_Man_t * pIfMan, If_Obj_t * pIfObj, Vec_Int_t * vCover );
extern Abc_Obj_t * Abc_NodeFromIfDec_rec( Abc_Ntk_t * pNtk, Abc_Ntk_t * pNtkNew, If_Man_t * pIfMan, If_Obj_t * pIfObj, Vec_Int_t * vCover );
static Hop_Obj_t * Abc_NodeIfToHop( Hop_Man_t * pHopMan, If_Man_t * pIfMan, If_Obj_t * pIfObj );
static Vec_Ptr_t * Abc_NtkFindGoodOrder( Abc_Ntk_t * pNtk );

extern void Abc_NtkBddReorder( Abc_Ntk_t * pNtk, int fVerbose );
extern void Abc_NtkBidecResyn( Abc_Ntk_t * pNtk, int fVerbose );
 
////////////////////////////////////////////////////////////////////////
///                     FUNCTION DEFINITIONS                         ///
////////////////////////////////////////////////////////////////////////

/**Function*************************************************************

  Synopsis    [Interface with the FPGA mapping package.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void If_ManComputeSwitching( If_Man_t * pIfMan )
{
    abctime clk = Abc_Clock();
    Gia_Man_t * pNew;
    Vec_Int_t * vCopy;
    If_Obj_t * pIfObj;
    int i;
    assert( pIfMan->vSwitching == NULL );
    // create the new manager
    pNew = Gia_ManStart( If_ManObjNum(pIfMan) );
    vCopy = Vec_IntAlloc( If_ManObjNum(pIfMan) );
    // constant and inputs
    Vec_IntPush( vCopy, 1 );
    If_ManForEachCi( pIfMan, pIfObj, i )
        Vec_IntPush( vCopy, Gia_ManAppendCi(pNew) );
    // internal nodes
    If_ManForEachNode( pIfMan, pIfObj, i )
    {
        int iLit0 = Abc_LitNotCond( Vec_IntEntry(vCopy, If_ObjFanin0(pIfObj)->Id), If_ObjFaninC0(pIfObj) );
        int iLit1 = Abc_LitNotCond( Vec_IntEntry(vCopy, If_ObjFanin1(pIfObj)->Id), If_ObjFaninC1(pIfObj) );
        Vec_IntPush( vCopy, Gia_ManAppendAnd(pNew, iLit0, iLit1) );
    }
    // outputs
    If_ManForEachCo( pIfMan, pIfObj, i )
    {
        int iLit0 = Abc_LitNotCond( Vec_IntEntry(vCopy, If_ObjFanin0(pIfObj)->Id), If_ObjFaninC0(pIfObj) );
        Vec_IntPush( vCopy, Gia_ManAppendCo(pNew, iLit0) );
    }
    assert( Vec_IntSize(vCopy) == If_ManObjNum(pIfMan) );
    Vec_IntFree( vCopy );
    // compute switching activity
    pIfMan->vSwitching = Gia_ManComputeSwitchProbs( pNew, 48, 16, 0 );
    Gia_ManStop( pNew );
    if ( pIfMan->pPars->fVerbose )
        Abc_PrintTime( 1, "Computing switching activity", Abc_Clock() - clk );
}

int Abc_IfNtkFindMaxLevel(Abc_Ntk_t *pNtk)
{
    int maxLevel = 0;
    Abc_Obj_t *pNode;
    int i;
    // find maxLevel
    Abc_NtkForEachNode(pNtk, pNode, i)
    {
        if (Abc_ObjLevel(pNode) > maxLevel)
        {
            maxLevel = Abc_ObjLevel(pNode);
        }
    }
    return maxLevel;
}

Vec_Ptr_t * Abc_IfGetNodesAtLevel( Abc_Ntk_t *pNtk, u_int32_t level )
{
    Vec_Ptr_t * nodesAtLevel = Vec_PtrAlloc(0);
    Abc_Obj_t * pNode;
    int i;

    Abc_NtkForEachCi(pNtk, pNode, i)
    {
        if( Abc_ObjLevel( pNode ) == level )
        {
            Vec_PtrPush(nodesAtLevel, pNode);
        }
    }
    Abc_NtkForEachNode(pNtk, pNode, i)
    {
        if( Abc_ObjLevel( pNode ) == level )
        {
            Vec_PtrPush(nodesAtLevel, pNode);
        }
    }
    return nodesAtLevel;
}

typedef struct
{
    Abc_Obj_t * source;
    Abc_Obj_t * target;
} Edge;

typedef struct
{
    double x, y;
} NodePos;

NodePos Abc_IfGetNodePos( Abc_Ntk_t * pNtk , Abc_Obj_t * n)
{
    NodePos result;

    result.x = n->iTemp;
    result.y = n->Level;

    return result;
}

int Abc_IfIsStraighLineCrossing( Abc_Ntk_t * pNtk, Abc_Obj_t * src1, Abc_Obj_t * tgt1, Abc_Obj_t * src2, Abc_Obj_t * tgt2 )
{
    NodePos p_src1 = Abc_IfGetNodePos( pNtk, src1 );
    NodePos p_tgt1 = Abc_IfGetNodePos( pNtk, tgt1 );
    NodePos p_src2 = Abc_IfGetNodePos( pNtk, src2 );
    NodePos p_tgt2 = Abc_IfGetNodePos( pNtk, tgt2 );

    double s1_x = p_tgt1.x - p_src1.x;
    double s1_y = p_tgt1.y - p_src1.y;
    double s2_x = p_tgt2.x - p_src2.x;
    double s2_y = p_tgt2.y - p_src2.y;

    double s = (-s1_y * (p_src1.x - p_src2.x) + s1_x * (p_src1.y - p_src2.y)) / (-s2_x * s1_y + s1_x * s2_y);
    double t = (s2_x * (p_src1.y - p_src2.y) - s2_y * (p_src1.x - p_src2.x)) / (-s2_x * s1_y + s1_x * s2_y);

    // do not use >= / <= here to not consider lines that share the same endpoints as a crossing
    return (s > 0 && s < 1 && t > 0 && t < 1);
}

size_t Abc_NumCrossingsBetweenLevels( Abc_Ntk_t * pNtk, u_int32_t l1, u_int32_t l2 )
{
    // There are no crossings with itself
    if (l1 == l2) {
        return 0;
    }
    // Make sure that r1 is the lower level
    if (l1 > l2) {
        u_int32_t temp = l1;
        l1 = l2;
        l2 = temp;
    }

    assert(l1 + 1 == l2 && "r1 and r2 must be consecutive levels");

    size_t crossings = 0;

    // Gather all edges between levels r1 and r2
    Vec_Ptr_t* edges = Vec_PtrAlloc(0);
    Abc_Obj_t *pNode, *pFanin;
    int i, j;

    Vec_Ptr_t* nodesAtLevel2 =  Abc_IfGetNodesAtLevel( pNtk, l2 ) ;

    Vec_PtrForEachEntry(Abc_Obj_t *, nodesAtLevel2, pNode, i) {
        Abc_ObjForEachFanin(pNode, pFanin, j) {
            if(Abc_ObjLevel(pFanin) == l1) {
                Edge* edge = (Edge*)malloc(sizeof(Edge));
                edge->source = pFanin;
                edge->target = pNode;
                Vec_PtrPush(edges, edge);

                /*printf("Source: %i, Target: %i\n", pFanin->Id, pNode->Id);
                printf("Edge level: %d, Source (x=%d, y=%d), Target (x=%d, y=%d)\n",
                       l2,
                       edge->source->iTemp,
                       Abc_ObjLevel(edge->source),
                       edge->target->iTemp,
                       Abc_ObjLevel(edge->target));*/
            }
        }
    }

    int nEdges = Vec_PtrSize(edges);

    if(nEdges >= 2) {
        for(int i = 0; i < nEdges - 1; i++) {
            for(int j = i + 1; j < nEdges; j++) {
                Edge* edge1 = (Edge*)Vec_PtrEntry(edges, i);
                Edge* edge2 = (Edge*)Vec_PtrEntry(edges, j);
                if(Abc_IfIsStraighLineCrossing( pNtk, edge1->source, edge1->target, edge2->source, edge2->target ) )
                {
                    crossings++;
                }
            }
        }
    }

    // free memory
    for(int k = 0; k < nEdges; k++) {
        free(Vec_PtrEntry(edges, k));
    }
    Vec_PtrFree( edges );
    Vec_PtrFree( nodesAtLevel2 );

    return crossings;
}

size_t Abc_IfNumCrossingsWithAdjacentLevels( Abc_Ntk_t * pNtk, u_int32_t l )
{
    size_t crossings = 0;
    int maxLevel = Abc_IfNtkFindMaxLevel( pNtk );
    // If r is not the last level, add crossings with level r + 1
    if (l != maxLevel) {
        crossings += Abc_NumCrossingsBetweenLevels(pNtk, l, l + 1);
    }

    // If r is not 0, there are no crossings with the previous level
    if (l != 0) {
        // Add crossings with level r - 1
        crossings += Abc_NumCrossingsBetweenLevels(pNtk, l - 1, l);
    }

    return crossings;
}

void InsertBufferNodes( Abc_Ntk_t * pNtk )
{
    Abc_Obj_t *pNode;
    int i;
    // Abc_Ntk_t * pNtkNew = Abc_NtkDup( pNtk );

    // Set levels for all nodes
    Abc_NtkLevel( pNtk );

    Abc_NtkForEachNode( pNtk, pNode, i )
    {
        Abc_Obj_t *pFanin;
        int j;
        Abc_ObjForEachFanin( pNode, pFanin, j )
        {
            for ( int lv = Abc_ObjLevel(pFanin) + 1; lv < Abc_ObjLevel(pNode); lv++ )
            {
                // create new buffer node
                Abc_Obj_t *pBufferNode = Abc_NtkCreateNode( pNtk );

                // Disconnect the original edge
                Abc_ObjPatchFanin( pNode, pFanin, pBufferNode );

                // Connect the dummy buffer to the fanin
                Abc_ObjAddFanin( pBufferNode, pFanin );

                // update pFanin of the Node
                pFanin = Abc_ObjFanin( pNode, j );
            }
        }
    }
    // Set levels again for all nodes
    pNtk->LevelMax = Abc_NtkLevel( pNtk );
    if ( !Abc_NtkCheck( pNtk ) )
    {
        printf("Error: Buffered Ntk doesnt CHECK\n");
    }
}

void Abc_IfComputeRanks(Abc_Ntk_t* pNtk, int** rankArray)
{
    Abc_Obj_t* pNode;
    Abc_Obj_t* pObjDfs;
    // Vec_Ptr_t * vNodesDfs; Hallo
    int i, iDfs, currLevel;
    int maxLevel = Abc_IfNtkFindMaxLevel( pNtk );

    // create and initialize rankArray
    *rankArray = (int*)calloc(maxLevel + 1, sizeof(int));

    // compute ranks
    // vNodesDfs = Abc_NtkDfs(pNtk, 0);
    Abc_NtkForEachCi(pNtk, pNode, i)
    {
        pNode->iTemp = pNode->Id - 1 ;
        /*printf("NodeID: %i \n", pNode->Id);
        printf("Level: %i \n", pNode->Level);*/

    }
    Abc_NtkForEachNode( pNtk, pObjDfs, iDfs )
    {
        // get the current level currLevel
        currLevel = Abc_ObjLevel( pObjDfs );
        // rankArray[currLevel] holds the next available rank on level currLevel
        // use iTemp to temporarily store the rank of the node
        pObjDfs->iTemp = (*rankArray)[currLevel]++;
        /*printf("NodeID: %i \n", pObjDfs->Id);
        printf("Level: %i \n", pObjDfs->Level);*/
    }
    // Vec_PtrFree( vNodesDfs );
}

size_t Abc_IfCount_Crossings( Abc_Ntk_t * pNtk )
{
    size_t crossings = 0;
    int maxLevel = Abc_IfNtkFindMaxLevel( pNtk );
    u_int32_t l;
    for(l = 0; l <= maxLevel; l += 2) {
        crossings += Abc_IfNumCrossingsWithAdjacentLevels(pNtk, l);
    }
    // printf("CrossingNum: %zu\n", crossings);
    return crossings;
}

size_t Abc_IfComputeCrossingNum( Abc_Ntk_t * pNtk )
{
    //Abc_Ntk_t * pNtkNew;

    // insert buffers
    InsertBufferNodes( pNtk );
    if( !pNtk )
    {
        printf("Error while inserting Buffers\n");
    }

    // create and initialize rankArray
    int * rankArray;
    Abc_IfComputeRanks( pNtk, &rankArray );

    // compute crossings
    size_t crossings;

    // ntk.depth() in mockturtle is equivalent to maxLevel in ABC
    crossings = Abc_IfCount_Crossings( pNtk );

    free(rankArray); // remember to free the allocated memory after you're done

    //Abc_NtkDelete( pNtkNew );

    return crossings;
}

void If_ManComputeCrossings( Abc_Ntk_t * pNtk )
{
    Abc_Ntk_t * pNtkOpt = Abc_NtkDup( pNtk );

    size_t crossings_b = Abc_IfComputeCrossingNum( pNtkOpt );
    printf("CrossingNum: %zu\n", crossings_b);

    // Compute the ranks so that they are on the average positions between nodes


   /* size_t crossings_a = Abc_IfComputeCrossingNum( pNtkOpt );
    printf("CrossingRed: %zu\n", crossings_a-crossings_b);*/

    Abc_NtkDelete( pNtkOpt );
}

Abc_Ntk_t * If_ManReduceCrossings( If_Man_t * pIfMan, Abc_Ntk_t * pNtk )
{
    ProgressBar * pProgress;
    Abc_Ntk_t * pNtkNew;
    Abc_Obj_t * pNode, * pNodeNew;
    Vec_Int_t * vCover;
    int i, nDupGates;
    // create the new network
    if ( pIfMan->pPars->fUseBdds || pIfMan->pPars->fUseCnfs || pIfMan->pPars->fUseMv )
        pNtkNew = Abc_NtkStartFrom( pNtk, ABC_NTK_LOGIC, ABC_FUNC_BDD );
    else if ( pIfMan->pPars->fUseSops || pIfMan->pPars->fUserSesLib || pIfMan->pPars->nGateSize > 0 )
        pNtkNew = Abc_NtkStartFrom( pNtk, ABC_NTK_LOGIC, ABC_FUNC_SOP );
    else
        pNtkNew = Abc_NtkStartFrom( pNtk, ABC_NTK_LOGIC, ABC_FUNC_AIG );
    // prepare the mapping manager
    If_ManCleanNodeCopy( pIfMan );
    If_ManCleanCutData( pIfMan );
    // make the mapper point to the new network
    If_ObjSetCopy( If_ManConst1(pIfMan), Abc_NtkCreateNodeConst1(pNtkNew) );
    Abc_NtkForEachCi( pNtk, pNode, i )
        If_ObjSetCopy( If_ManCi(pIfMan, i), pNode->pCopy );

    // process the nodes in topological order
    vCover = Vec_IntAlloc( 1 << 16 );
    pProgress = Extra_ProgressBarStart( stdout, Abc_NtkCoNum(pNtk) );
    Abc_NtkForEachCo( pNtk, pNode, i )
    {
        Extra_ProgressBarUpdate( pProgress, i, "Final" );
        pNodeNew = Abc_NodeFromIfDec_rec( pNtk, pNtkNew, pIfMan, If_ObjFanin0(If_ManCo(pIfMan, i)), vCover );
        pNodeNew = Abc_ObjNotCond( pNodeNew, If_ObjFaninC0(If_ManCo(pIfMan, i)) );
        Abc_ObjAddFanin( pNode->pCopy, pNodeNew );
    }
    Extra_ProgressBarStop( pProgress );
    Vec_IntFree( vCover );

    // remove the constant node if not used
    pNodeNew = (Abc_Obj_t *)If_ObjCopy( If_ManConst1(pIfMan) );
    if ( Abc_ObjFanoutNum(pNodeNew) == 0 && !Abc_ObjIsNone(pNodeNew) )
        Abc_NtkDeleteObj( pNodeNew );
    // minimize the node
    if ( pIfMan->pPars->fUseBdds || pIfMan->pPars->fUseCnfs || pIfMan->pPars->fUseMv )
        Abc_NtkSweep( pNtkNew, 0 );
    if ( pIfMan->pPars->fUseBdds )
        Abc_NtkBddReorder( pNtkNew, 0 );
    // decouple the PO driver nodes to reduce the number of levels
    nDupGates = Abc_NtkLogicMakeSimpleCos( pNtkNew, !pIfMan->pPars->fUseBuffs );
    if ( nDupGates && pIfMan->pPars->fVerbose && !Abc_FrameReadFlag("silentmode") )
    {
        if ( pIfMan->pPars->fUseBuffs )
            printf( "Added %d buffers/inverters to decouple the CO drivers.\n", nDupGates );
        else
            printf( "Duplicated %d gates to decouple the CO drivers.\n", nDupGates );
    }
    return pNtkNew;
}
/**Function*************************************************************

  Synopsis    []

  Description [Traversal functions to find shared nodes in AIG.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
/*Helper function*/
void If_ManCleanMarkS( If_Man_t * p )
{
    If_Obj_t * pObj;
    int i;
    If_ManForEachObj( p, pObj, i )
        pObj->fSpec = 0;
}
int CountSharedLeaves( If_Cut_t * pC0, If_Cut_t * pC1 )
{
    int nSzC0 = pC0->nLeaves;
    int nSzC1 = pC1->nLeaves;
    int sharedLeaves = 0;
    for ( int l0 = 0; l0 < nSzC0; l0++ )
    {
        for ( int l1 = 0; l1 < nSzC1; l1++ )
        {
            if ( pC0->pLeaves[l0] == pC1->pLeaves[l1] )
            {
                sharedLeaves++;
                break;
            }
        }
    }
    return sharedLeaves;
}
/*void Decompose( Abc_Ntk_t * pNtk, If_Man_t * pIfMan, If_Cut_t * pC0, If_Cut_t * pC1 )
{

}
int Abc_TraverseandCompareNodesIf_rec( Abc_Ntk_t * pNtk, If_Man_t * pIfMan, If_Obj_t * pIfObj, If_Cut_t * pC0 )
{
    If_Cut_t * pC1;
    If_Obj_t * pIfLeaf;
    int i;
    if( pIfObj->Type != IF_AND || pIfObj->fVisit == 1 || pIfObj->fSpec )
        return 1;

    pC1 = If_ObjCutBest( pIfObj );

    If_CutForEachLeaf(pIfMan, pC1, pIfLeaf, i)
    {
        Abc_TraverseandCompareNodesIf_rec( pNtk, pIfMan, pIfLeaf, pC0 ); // recursive call
    }

    printf("Second Id: %i\n", pIfObj->Id);
    int sharedLeaves = CountSharedLeaves( pC0, pC1 );
    if( sharedLeaves > 1 )
        printf("Number of shared leaves: %d\n", sharedLeaves);
    if( sharedLeaves == 4 )
        Decompose( pNtk, pIfMan, pC0, pC1 );

    pIfObj->fSpec = 1;

    return 1;
}

void compareCuts( Abc_Ntk_t * pNtk, If_Man_t * pIfMan, If_Cut_t * pC )
{
    Abc_Obj_t * pNode;
    int i;
    If_ManCleanMarkS( pIfMan );
    Abc_NtkForEachCo( pNtk, pNode, i )
    {
        Abc_TraverseandCompareNodesIf_rec( pNtk, pIfMan, If_ObjFanin0(If_ManCo(pIfMan, i)), pC );
    }
}*/

int Abc_TraverseNodesIfTwice_rec( Abc_Ntk_t * pNtk, If_Man_t * pIfMan, If_Obj_t * pIfObj )
{
    If_Cut_t * pC;
    If_Obj_t * pIfLeaf;
    int i;
    if( (pIfObj->Type != IF_AND) || pIfObj->fVisit == 1 )
        return 1;
    pC = If_ObjCutBest( pIfObj );

    If_CutForEachLeaf(pIfMan, pC, pIfLeaf, i)
    {
        Abc_TraverseNodesIfTwice_rec( pNtk, pIfMan, pIfLeaf ); // recursive call
    }
    printf("Id: %i\n", pIfObj->Id);
    pIfObj->fVisit = 1;
    // Last cut can be skipped since there is no cut to compare it to
    if ( pC == If_ObjCutBest( If_ObjFanin0( If_ManCo( pIfMan, Vec_PtrSize( pIfMan->vCos) - 1 ) ) ) )
    {
        return 0;
    }
    // compareCuts( pNtk, pIfMan, pC ); // compare the leaves

    return 0;
}
If_Obj_t * Abc_TraverseNodesIf_rec( Abc_Ntk_t * pNtk, If_Man_t * pIfMan, If_Obj_t * pIfObj, If_Cut_t * pC0 )
{
    If_Cut_t * pC1;
    If_Obj_t * pIfLeaf, * nodeResult;
    int i;

    if( (pIfObj->Type != IF_AND) || pIfObj->fVisit == 1 )
        return NULL;

    pC1 = If_ObjCutBest( pIfObj );

    If_CutForEachLeaf(pIfMan, pC1, pIfLeaf, i)
    {
        nodeResult = Abc_TraverseNodesIf_rec( pNtk, pIfMan, pIfLeaf, pC0 ); // recursive call
        if (nodeResult != NULL && pC0 != pC1) {
            return nodeResult; // break and return if a match is found
        }
    }

    int sharedLeaves = CountSharedLeaves( pC0, pC1 );
    if( sharedLeaves == 4 && pC1 != pC0 )
    {
        // printf("Number of shared leaves: %d\n", sharedLeaves);
        return pIfObj;
    }

    pIfObj->fVisit = 1;
    // Last cut can be skipped since there is no cut to compare it to
    if ( pC1 == If_ObjCutBest( If_ObjFanin0( If_ManCo( pIfMan, Vec_PtrSize( pIfMan->vCos) - 1 ) ) ) )
    {
        return NULL;
    }

    return NULL;
}
/**Function*************************************************************

  Synopsis    [Create a Mapping Manager with decomposed LUTs to decrease crossings.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void LutDecomposition( If_Man_t * pIfMan, Abc_Ntk_t * pNtk )
{
    If_Obj_t * pIfObj;
    Abc_Obj_t * pNode;
    int i;
    // go through pIfMan
    printf("Ein Aufruf\n");
    If_ManCleanMarkV( pIfMan );

    Abc_NtkForEachCo( pNtk, pNode, i )
    {
        printf("Iteration0 PO\n");
        Abc_TraverseNodesIfTwice_rec( pNtk, pIfMan, If_ObjFanin0(If_ManCo(pIfMan, i)) );
    }
    // find LUTs with the same shared inputs
    // find all decompositions of these LUTs
    // if it can find a valid solution accept it, if not try another LUT pair.
}

/**Function*************************************************************

  Synopsis    [Interface with the FPGA mapping package.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
Abc_Ntk_t * Abc_NtkIf( Abc_Ntk_t * pNtk, If_Par_t * pPars )
{
    Abc_Ntk_t * pNtkNew, * pTemp;
    If_Man_t * pIfMan;

    assert( Abc_NtkIsStrash(pNtk) );

    // get timing information
    pPars->pTimesArr = Abc_NtkGetCiArrivalFloats(pNtk);
    pPars->pTimesReq = Abc_NtkGetCoRequiredFloats(pNtk);

    // update timing info to reflect logic level
    if ( (pPars->fDelayOpt || pPars->fDsdBalance || pPars->fUserRecLib || pPars->fUserSesLib || pPars->fUserLutDec) && pNtk->pManTime )
    {
        int c;
        if ( pNtk->AndGateDelay == 0.0 )
        {
            if ( Abc_FrameReadLibGen() )
                pNtk->AndGateDelay = Mio_LibraryReadDelayAigNode((Mio_Library_t *)Abc_FrameReadLibGen());
            if ( pNtk->AndGateDelay == 0.0 )
            {
                pNtk->AndGateDelay = 1.0;
                printf( "The AIG-node delay is not set. Assuming unit-delay.\n" );
            }
        }
        for ( c = 0; c < Abc_NtkCiNum(pNtk); c++ )
            pPars->pTimesArr[c] /= pNtk->AndGateDelay;
        for ( c = 0; c < Abc_NtkCoNum(pNtk); c++ )
            pPars->pTimesReq[c] /= pNtk->AndGateDelay;
    }

    // set the latch paths
    if ( pPars->fLatchPaths && pPars->pTimesArr )
    {
        int c;
        for ( c = 0; c < Abc_NtkPiNum(pNtk); c++ )
            pPars->pTimesArr[c] = -ABC_INFINITY;
    }

    // create FPGA mapper
    pIfMan = Abc_NtkToIf( pNtk, pPars );    
    if ( pIfMan == NULL )
        return NULL;
    if ( pPars->fPower )
        If_ManComputeSwitching( pIfMan );

    // create DSD manager
    if ( pPars->fUseDsd )
    {
        If_DsdMan_t * p = (If_DsdMan_t *)Abc_FrameReadManDsd();
        assert( pPars->nLutSize <= If_DsdManVarNum(p) );
        assert( (pPars->pLutStruct == NULL && If_DsdManLutSize(p) == 0) || (pPars->pLutStruct && pPars->pLutStruct[0] - '0' == If_DsdManLutSize(p)) );
        pIfMan->pIfDsdMan = (If_DsdMan_t *)Abc_FrameReadManDsd();
        if ( pPars->fDsdBalance )
            If_DsdManAllocIsops( pIfMan->pIfDsdMan, pPars->nLutSize );
    }

    // perform FPGA mapping
    if ( !If_ManPerformMapping( pIfMan ) )
    {
        If_ManStop( pIfMan );
        return NULL;
    }

    // transform the result of mapping into the new network
    // pNtkNew = Abc_NtkFromIf( pIfMan, pNtk );
    pNtkNew = If_ManReduceCrossings( pIfMan, pNtk );
    if ( pNtkNew == NULL )
        return NULL;

    If_ManStop( pIfMan );
    if ( pPars->fDelayOpt || pPars->fDsdBalance || pPars->fUserRecLib )
    {
        pNtkNew = Abc_NtkStrash( pTemp = pNtkNew, 0, 0, 0 );
        Abc_NtkDelete( pTemp );
    }
    else if ( pPars->fBidec && pPars->nLutSize <= 8 )
        Abc_NtkBidecResyn( pNtkNew, 0 );

    // duplicate EXDC
    if ( pNtk->pExdc )
        pNtkNew->pExdc = Abc_NtkDup( pNtk->pExdc );
    // make sure that everything is okay
    if ( !Abc_NtkCheck( pNtkNew ) )
    {
        printf( "Abc_NtkIf: The network check has failed.\n" );
        Abc_NtkDelete( pNtkNew );
        return NULL;
    }
    return pNtkNew;
}

/**Function*************************************************************

  Synopsis    [Load the network into FPGA manager.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
static inline If_Obj_t * Abc_ObjIfCopy( Abc_Obj_t * pNode ) { return (If_Obj_t *)pNode->pCopy; } 
If_Man_t * Abc_NtkToIf( Abc_Ntk_t * pNtk, If_Par_t * pPars )
{
    ProgressBar * pProgress;
    If_Man_t * pIfMan;
    Vec_Ptr_t * vNodes;
    Abc_Obj_t * pNode, * pPrev;
    int i;

    assert( Abc_NtkIsStrash(pNtk) );

    // start the mapping manager and set its parameters
    pIfMan = If_ManStart( pPars );
    pIfMan->pName = Abc_UtilStrsav( Abc_NtkName(pNtk) );

    // print warning about excessive memory usage
    if ( 1.0 * Abc_NtkObjNum(pNtk) * pIfMan->nObjBytes / (1<<30) > 1.0 )
        printf( "Warning: The mapper will allocate %.1f GB for to represent the subject graph with %d AIG nodes.\n", 
            1.0 * Abc_NtkObjNum(pNtk) * pIfMan->nObjBytes / (1<<30), Abc_NtkObjNum(pNtk) );

    // create PIs and remember them in the old nodes
    Abc_NtkCleanCopy( pNtk );
    Abc_AigConst1(pNtk)->pCopy = (Abc_Obj_t *)If_ManConst1( pIfMan );
    Abc_NtkForEachCi( pNtk, pNode, i )
    {
        If_Obj_t * pIfObj = If_ManCreateCi( pIfMan );
        pNode->pCopy = (Abc_Obj_t *)pIfObj;
        // transfer logic level information
        Abc_ObjIfCopy(pNode)->Level = pNode->Level;
        // mark the largest level
        if ( pIfMan->nLevelMax < (int)pIfObj->Level )
            pIfMan->nLevelMax = (int)pIfObj->Level;
    }

    // load the AIG into the mapper
    pProgress = Extra_ProgressBarStart( stdout, Abc_NtkObjNumMax(pNtk) );
    vNodes = Abc_AigDfs( pNtk, 0, 0 );
    Vec_PtrForEachEntry( Abc_Obj_t *, vNodes, pNode, i )
    {
        Extra_ProgressBarUpdate( pProgress, i, "Initial" );
        // add the node to the mapper
        pNode->pCopy = (Abc_Obj_t *)If_ManCreateAnd( pIfMan, 
            If_NotCond( Abc_ObjIfCopy(Abc_ObjFanin0(pNode)), Abc_ObjFaninC0(pNode) ), 
            If_NotCond( Abc_ObjIfCopy(Abc_ObjFanin1(pNode)), Abc_ObjFaninC1(pNode) ) );
        // set up the choice node
        if ( Abc_AigNodeIsChoice( pNode ) )
        {
            Abc_Obj_t * pEquiv;
//            int Counter = 0;
            assert( If_ObjId(Abc_ObjIfCopy(pNode)) > If_ObjId(Abc_ObjIfCopy(Abc_ObjEquiv(pNode))) );
            for ( pPrev = pNode, pEquiv = Abc_ObjEquiv(pPrev); pEquiv; pPrev = pEquiv, pEquiv = Abc_ObjEquiv(pPrev) )
                If_ObjSetChoice( Abc_ObjIfCopy(pPrev), Abc_ObjIfCopy(pEquiv) );//, Counter++;
//            printf( "%d ", Counter );
            If_ManCreateChoice( pIfMan, Abc_ObjIfCopy(pNode) );
        }
    }
    Extra_ProgressBarStop( pProgress );
    Vec_PtrFree( vNodes );

    // set the primary outputs without copying the phase
    Abc_NtkForEachCo( pNtk, pNode, i )
        pNode->pCopy = (Abc_Obj_t *)If_ManCreateCo( pIfMan, If_NotCond( Abc_ObjIfCopy(Abc_ObjFanin0(pNode)), Abc_ObjFaninC0(pNode) ) );
    return pIfMan;
}


/**Function*************************************************************

  Synopsis    [Box mapping procedures.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
static inline void Abc_MapBoxSetPrevNext( Vec_Ptr_t * vDrivers, Vec_Int_t * vMapIn, Vec_Int_t * vMapOut, int Id )
{
    Abc_Obj_t * pNode;
    pNode = (Abc_Obj_t *)Vec_PtrEntry(vDrivers, Id+2);
    Vec_IntWriteEntry( vMapIn, Abc_ObjId(Abc_ObjFanin0(pNode)), Id );
    pNode = (Abc_Obj_t *)Vec_PtrEntry(vDrivers, Id+4);
    Vec_IntWriteEntry( vMapOut, Abc_ObjId(Abc_ObjFanin0(pNode)), Id );
}
static inline int Abc_MapBox2Next( Vec_Ptr_t * vDrivers, Vec_Int_t * vMapIn, Vec_Int_t * vMapOut, int Id )
{
    Abc_Obj_t * pNode = (Abc_Obj_t *)Vec_PtrEntry(vDrivers, Id+4);
    return Vec_IntEntry( vMapIn, Abc_ObjId(Abc_ObjFanin0(pNode)) );
}
static inline int Abc_MapBox2Prev( Vec_Ptr_t * vDrivers, Vec_Int_t * vMapIn, Vec_Int_t * vMapOut, int Id )
{
    Abc_Obj_t * pNode = (Abc_Obj_t *)Vec_PtrEntry(vDrivers, Id+2);
    return Vec_IntEntry( vMapOut, Abc_ObjId(Abc_ObjFanin0(pNode)) );
}

/**Function*************************************************************

  Synopsis    [Creates the mapped network.]

  Description [Assuming the copy field of the mapped nodes are NULL.]
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
Abc_Ntk_t * Abc_NtkFromIf( If_Man_t * pIfMan, Abc_Ntk_t * pNtk )
{
    ProgressBar * pProgress;
    Abc_Ntk_t * pNtkNew;
    Abc_Obj_t * pNode, * pNodeNew;
    Vec_Int_t * vCover;
    int i, nDupGates;
    // create the new network
    if ( pIfMan->pPars->fUseBdds || pIfMan->pPars->fUseCnfs || pIfMan->pPars->fUseMv )
        pNtkNew = Abc_NtkStartFrom( pNtk, ABC_NTK_LOGIC, ABC_FUNC_BDD );
    else if ( pIfMan->pPars->fUseSops || pIfMan->pPars->fUserSesLib || pIfMan->pPars->nGateSize > 0 )
        pNtkNew = Abc_NtkStartFrom( pNtk, ABC_NTK_LOGIC, ABC_FUNC_SOP );
    else
        pNtkNew = Abc_NtkStartFrom( pNtk, ABC_NTK_LOGIC, ABC_FUNC_AIG );
    // prepare the mapping manager
    If_ManCleanNodeCopy( pIfMan );
    If_ManCleanCutData( pIfMan );
    // make the mapper point to the new network
    If_ObjSetCopy( If_ManConst1(pIfMan), Abc_NtkCreateNodeConst1(pNtkNew) );
    Abc_NtkForEachCi( pNtk, pNode, i )
        If_ObjSetCopy( If_ManCi(pIfMan, i), pNode->pCopy );

    // process the nodes in topological order
    vCover = Vec_IntAlloc( 1 << 16 );
    pProgress = Extra_ProgressBarStart( stdout, Abc_NtkCoNum(pNtk) );
    Abc_NtkForEachCo( pNtk, pNode, i )
    {
        Extra_ProgressBarUpdate( pProgress, i, "Final" );
        pNodeNew = Abc_NodeFromIf_rec( pNtkNew, pIfMan, If_ObjFanin0(If_ManCo(pIfMan, i)), vCover );
        pNodeNew = Abc_ObjNotCond( pNodeNew, If_ObjFaninC0(If_ManCo(pIfMan, i)) );
        Abc_ObjAddFanin( pNode->pCopy, pNodeNew );
    }
    Extra_ProgressBarStop( pProgress );
    Vec_IntFree( vCover );

    // remove the constant node if not used
    pNodeNew = (Abc_Obj_t *)If_ObjCopy( If_ManConst1(pIfMan) );
    if ( Abc_ObjFanoutNum(pNodeNew) == 0 && !Abc_ObjIsNone(pNodeNew) )
        Abc_NtkDeleteObj( pNodeNew );
    // minimize the node
    if ( pIfMan->pPars->fUseBdds || pIfMan->pPars->fUseCnfs || pIfMan->pPars->fUseMv )
        Abc_NtkSweep( pNtkNew, 0 );
    if ( pIfMan->pPars->fUseBdds )
        Abc_NtkBddReorder( pNtkNew, 0 );
    // decouple the PO driver nodes to reduce the number of levels
    nDupGates = Abc_NtkLogicMakeSimpleCos( pNtkNew, !pIfMan->pPars->fUseBuffs );
    if ( nDupGates && pIfMan->pPars->fVerbose && !Abc_FrameReadFlag("silentmode") )
    {
        if ( pIfMan->pPars->fUseBuffs )
            printf( "Added %d buffers/inverters to decouple the CO drivers.\n", nDupGates );
        else
            printf( "Duplicated %d gates to decouple the CO drivers.\n", nDupGates );
    }
    return pNtkNew;
}

/**Function*************************************************************

  Synopsis    [Rebuilds GIA from mini AIG.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
Hop_Obj_t * Abc_NodeBuildFromMiniInt( Hop_Man_t * pMan, Vec_Int_t * vAig, int nLeaves )
{
    assert( Vec_IntSize(vAig) > 0 );
    assert( Vec_IntEntryLast(vAig) < 2 );
    if ( Vec_IntSize(vAig) == 1 ) // const
    {
        assert( nLeaves == 0 );
        return Hop_NotCond( Hop_ManConst0(pMan), Vec_IntEntry(vAig, 0) );
    }
    if ( Vec_IntSize(vAig) == 2 ) // variable
    {
        assert( Vec_IntEntry(vAig, 0) == 0 );
        assert( nLeaves == 1 );
        return Hop_NotCond( Hop_IthVar(pMan, 0), Vec_IntEntry(vAig, 1) );
    }
    else
    {
        int i, iVar0, iVar1, iLit0, iLit1;
        Hop_Obj_t * piLit0, * piLit1, * piLit = NULL;
        assert( Vec_IntSize(vAig) & 1 );
        Vec_IntForEachEntryDouble( vAig, iLit0, iLit1, i )
        {
            iVar0 = Abc_Lit2Var( iLit0 );
            iVar1 = Abc_Lit2Var( iLit1 );
            piLit0 = Hop_NotCond( iVar0 < nLeaves ? Hop_IthVar(pMan, iVar0) : (Hop_Obj_t *)Vec_PtrEntry((Vec_Ptr_t *)vAig, iVar0 - nLeaves), Abc_LitIsCompl(iLit0) );
            piLit1 = Hop_NotCond( iVar1 < nLeaves ? Hop_IthVar(pMan, iVar1) : (Hop_Obj_t *)Vec_PtrEntry((Vec_Ptr_t *)vAig, iVar1 - nLeaves), Abc_LitIsCompl(iLit1) );
            piLit  = Hop_And( pMan, piLit0, piLit1 );
            assert( (i & 1) == 0 );
            Vec_PtrWriteEntry( (Vec_Ptr_t *)vAig, Abc_Lit2Var(i), piLit );  // overwriting entries
        }
        assert( i == Vec_IntSize(vAig) - 1 );
        piLit = Hop_NotCond( piLit, Vec_IntEntry(vAig, i) );
        Vec_IntClear( vAig ); // useless
        return piLit;
    }
}
Hop_Obj_t * Abc_NodeBuildFromMini( Hop_Man_t * pMan, If_Man_t * p, If_Cut_t * pCut, int fUseDsd )
{
    int Delay;
    if ( fUseDsd )
        Delay = If_CutDsdBalanceEval( p, pCut, p->vArray );
    else
        Delay = If_CutSopBalanceEval( p, pCut, p->vArray );
    assert( Delay >= 0 );
    return Abc_NodeBuildFromMiniInt( pMan, p->vArray, If_CutLeaveNum(pCut) );
}

/**Function*************************************************************

  Synopsis    [Implements decomposed LUT-structure of the cut.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void printBits(word input) {
    int bits_in_word = 8 * sizeof(word);
    for (int bit = 0; bit < bits_in_word; ++bit) {
        printf("%d", (input >> (bits_in_word - bit - 1)) & 1);
    }
    printf("\n");
}
void printBits2(word input) {
    int bits_in_word = 4;
    for (int bit = 0; bit < bits_in_word; ++bit) {
        printf("%d", (input >> (bits_in_word - bit - 1)) & 1);
    }
    printf("\n");
}
// Function to convert 16, 8, or 4-bit sequence into 64-bit sequence by repeating it
word convertTo64Bit(word smallTruth, int bits) {
    word result = 0;
    int repeats = 64 / bits;

    for (int i = 0; i < repeats; i++)
        result = (result << bits) | smallTruth;

    return result;
}

void swapAdjacentIterations(word *pTruth, int nWords, int iVarIndices[], int numIndices) {
    for (int i = 0; i < numIndices; i++) {
        Abc_TtSwapAdjacent(pTruth, nWords, iVarIndices[i]);
    }
}

void Abc_TtShift(word *pTruth, int nWords, int nVars, int nShifts) {
    if (nVars <= 1 || nShifts == 0) {
        return;
    }

    switch (nShifts) {
        case 1: {
            int indices[] = {2};
            swapAdjacentIterations(pTruth, nWords, indices, sizeof(indices) / sizeof(indices[0]));
            break;
        }
        case 2: {
            int indices[] = {1, 2, 1};
            swapAdjacentIterations(pTruth, nWords, indices, sizeof(indices) / sizeof(indices[0]));
            break;
        }
        case 3: {
            int indices[] = {0, 1, 2, 1, 0};
            swapAdjacentIterations(pTruth, nWords, indices, sizeof(indices) / sizeof(indices[0]));
            break;
        }

        default:
            printf("Invalid bit value\n");
            break;
    }
}
void Abc_TtCircularShift(word * pTruth, int nWords, int nVars)
{
    if(nVars <= 1) // If there is only 1 variable or none, simply return
    {
        return;
    }

    // For every variable from 0 to nVars-1, swap iVar with iVar+1
    for(int iVar=0; iVar < nVars-1; iVar++)
    {
        Abc_TtSwapAdjacent(pTruth, nWords, iVar);
    }
}
void Abc_TtLValue(word * pTruth, int nWords, int nVars)
{
    if(nVars <= 1) // If there is only 1 variable or none, simply return
    {
        return;
    }

    // For every variable from nVars-1 to 1, swap iVar with iVar-1
    for(int iVar = nVars-1; iVar > 0; iVar--)
    {
        Abc_TtSwapAdjacent(pTruth, nWords, iVar-1);
    }
}
int apply_operations( word lut4inp, word lut3inp, word lut2inp )
{
    word bit1 = 0xFF00 ;
    bit1 = convertTo64Bit(bit1, 16);
    word output = 0;
    int bits_in_word = 4;
    for (int bit = 0; bit < bits_in_word; ++bit) {
        int do_operation = (lut2inp >> (bits_in_word - bit - 1)) & 1;

        if ( do_operation )
        {
            switch (bit)
            {
                case 0:
                    output = output | ( bit1 & lut3inp );
                    break;

                case 1:
                    output = output | ( bit1 & ~lut3inp );
                    break;

                case 2:
                    output = output | ( ~bit1 & lut3inp );
                    break;

                case 3:
                    output = output | ( ~bit1 & ~lut3inp );
                    break;

                default:
                    printf("Invalid bit value\n");
                    break;
            }
        }
    }
    if( lut4inp == output )
    {
        return 1;
    }


    return 0;
}

// Function to generate the binary sequence and return it as an array
word* generate_binary_sequence(int num_bits, int* size) {
    *size = (1 << num_bits);
    word* sequence = (word*)malloc(*size * sizeof(word));
    sequence[0] = 0; // First entry is all 0s
    int idx = 1;
    for (int num_bits_set = 1; num_bits_set <= num_bits; ++num_bits_set) {
        for (int i = 0; i < (1 << num_bits); ++i) {
            int set_bits_count = 0;
            for (int j = 0; j < num_bits; ++j) {
                if ((i >> j) & 1) {
                    set_bits_count++;
                }
            }
            if (set_bits_count == num_bits_set) {
                sequence[idx++] = i;
            }
        }
    }
    return sequence;
}

int dec1( word four_inp_lut1, word four_inp_lut2 )
{
    int num_bits = 8;
    int size;
    word* sequence = generate_binary_sequence(num_bits, &size);
    word cur3inpLut, cur3inpLuttest;
    // find a three input LUT
    /*word* sequencetest = generate_binary_sequence(num_bits, &size);
    for (int i = 8; i < size; ++i)
    {
        cur3inpLuttest = convertTo64Bit(sequencetest[i], num_bits );
        printf("Before:\n");
        printBits(cur3inpLuttest);
        Abc_TtCircularShift( &cur3inpLuttest, 1, 4 );
        printf("After:\n");
        printBits(cur3inpLuttest);
        break;
    }*/
    for (int i = 0; i < size; ++i)
    {
        int found_dec = 0;
        cur3inpLut = convertTo64Bit(sequence[i], num_bits );
        // check for correct decomposition for the first four input LUT
        for (int j = 0; j < 16; j++)
        {
            found_dec = apply_operations(four_inp_lut1, cur3inpLut, (word) j);
            if ( found_dec )
            {
                /*printf("First Dec found: \n");
                printBits(cur3inpLut);*/
                // check for correct decomposition for the second four input LUT
                for (int k = 0; k < 16; k++)
                {
                    // Abc_TtCircularShift( &cur3inpLut, 1, 4 );
                    found_dec = apply_operations(four_inp_lut2, cur3inpLut, (word) k);
                    if ( found_dec )
                    {
                        printf("Four input LUTs: \n");
                        printBits(four_inp_lut1);
                        printBits(four_inp_lut2);
                        printf("Three input LUT: \n");
                        printBits(cur3inpLut);
                        printf("Two input LUTs: \n");
                        printBits((word) j);
                        printBits((word) k);

                        free(sequence);
                        return 1;
                    }
                }
            }

        }
    }

    // Free dynamically allocated memory
    free(sequence);
    return 0;
}


Hop_Obj_t * Abc_DecRecordToHop( Hop_Man_t * pMan, If_Man_t * pIfMan, If_Cut_t * pCutBest, If_Obj_t * pIfObj, Vec_Int_t * vCover )
{
    // get the truth table
    // perform LUT-decomposition and return the LUT-structure
    // convert the LUT-structure into a set of logic nodes in Abc_Ntk_t 

    // this is a placeholder, which takes the truth table and converts it into an AIG without LUT-decomposition
    extern Hop_Obj_t * Kit_TruthToHop( Hop_Man_t * pMan, unsigned * pTruth, int nVars, Vec_Int_t * vMemory );
    word * pTruth = If_CutTruthW(pIfMan, pCutBest);
    assert( !pIfMan->pPars->fUseTtPerm );
    return Kit_TruthToHop( (Hop_Man_t *)pMan, (unsigned *)pTruth, If_CutLeaveNum(pCutBest), vCover );
}
Hop_Obj_t * Abc_DecRecordToHopCros( Abc_Ntk_t * pNtk, Hop_Man_t * pMan, If_Man_t * pIfMan, If_Cut_t * pCutBest, If_Obj_t * pIfObj, Vec_Int_t * vCover )
{
    // get another truth table with shared leafs
    Abc_Obj_t * pNode;
    If_Obj_t * foundNode;
    If_Cut_t * pCutBest1;
    int i;
    printf("Function call\n");
    // go through pIfMan
    If_ManCleanMarkV( pIfMan );
    Abc_NtkForEachCo( pNtk, pNode, i )
    {
        foundNode = Abc_TraverseNodesIf_rec( pNtk, pIfMan, If_ObjFanin0(If_ManCo(pIfMan, i)), pCutBest );
        if (foundNode != NULL)
        {
            // printf("Ended search as node with sharedLeaves == 4 is found\n");
            break;
        }
    }
    pCutBest1 = If_ObjCutBest( foundNode );
    If_CutRotatePins( pIfMan, pCutBest1 );
    // depending on the order of leafs the tt has to be adjusted

    // int nSzC0 = pCutBest1->nLeaves;
    // int nSzC1 = pC1->nLeaves;
    /*for ( int l0 = 0; l0 < nSzC0; l0++ )
    {
        printf("Leaves: %i \n", pCutBest1->pLeaves[l0]);
    }*/
    // printf( "The nodes %i and %i have shared Leaves\n", pIfObj->Id, foundNode->Id );

    // get the truth table
    // perform LUT-decomposition and return the LUT-structure
    // convert the LUT-structure into a set of logic nodes in Abc_Ntk_t

    // this is a placeholder, which takes the truth table and converts it into an AIG without LUT-decomposition
    extern Hop_Obj_t * Kit_TruthToHop( Hop_Man_t * pMan, unsigned * pTruth, int nVars, Vec_Int_t * vMemory );
    word * pTruth1 = If_CutTruthW(pIfMan, pCutBest1);
    word pTruth1_copy = *pTruth1;

    word * pTruth0 = If_CutTruthW(pIfMan, pCutBest);
    word pTruth0_copy = *pTruth0;
    // both truth tables need to be in their non-complemented form
    if(  If_CutTruthIsCompl( pCutBest ) )
    {
        pTruth0_copy = ~pTruth0_copy;
    }
    if(  If_CutTruthIsCompl( pCutBest1 ) )
    {
        pTruth1_copy = ~ pTruth1_copy;
    }
    // word dcba
    // iteration 0 has cba as three input LUT, 1 = dba, 2 = dca, 3 = dcb.
    for(int iter_var=0; iter_var < 4; iter_var++)
    {
        // Apply the shift-to-front operation
        Abc_TtShift( &pTruth0_copy, 1, 4, iter_var);
        Abc_TtShift( &pTruth1_copy, 1, 4, iter_var);

        // Then call dec1 with the updated truth tables
        if( dec1(pTruth0_copy, pTruth1_copy) )
        {
            printf( "decomposition found: %i \n", iter_var );
            break;
        }
    }

    assert( !pIfMan->pPars->fUseTtPerm );
    return Kit_TruthToHop( (Hop_Man_t *)pMan, (unsigned *)pTruth0, If_CutLeaveNum(pCutBest), vCover );
}

Abc_Obj_t * Abc_NodeFromIfDec_rec( Abc_Ntk_t * pNtk, Abc_Ntk_t * pNtkNew, If_Man_t * pIfMan, If_Obj_t * pIfObj, Vec_Int_t * vCover )
{
    Abc_Obj_t * pNode, * pNodeNew;
    If_Cut_t * pCutBest;
    If_Obj_t * pIfLeaf;
    int i;
    // return if the result is known
    pNodeNew = (Abc_Obj_t *)If_ObjCopy( pIfObj );
    if ( pNodeNew )
        return pNodeNew;
    assert( pIfObj->Type == IF_AND );

    // get the parameters of the best cut
    pCutBest = If_ObjCutBest( pIfObj );
    if ( pIfMan->pPars->fUserSesLib )
    {
        // create the subgraph composed of Abc_Obj_t nodes based on the given cut
        Abc_Obj_t * pFanins[IF_MAX_FUNC_LUTSIZE];
        If_CutForEachLeaf( pIfMan, pCutBest, pIfLeaf, i )
            pFanins[i] = Abc_NodeFromIfDec_rec(pNtk, pNtkNew, pIfMan, pIfLeaf, vCover);
        pNodeNew = Abc_ExactBuildNode( If_CutTruthW(pIfMan, pCutBest), If_CutLeaveNum(pCutBest), If_CutArrTimeProfile(pIfMan, pCutBest), pFanins, pNtkNew );
        If_ObjSetCopy( pIfObj, pNodeNew );
        return pNodeNew;
    }
    // create a new node
    pNodeNew = Abc_NtkCreateNode( pNtkNew );

//    if ( pIfMan->pPars->pLutLib && pIfMan->pPars->pLutLib->fVarPinDelays )
    if ( !pIfMan->pPars->fDelayOpt && !pIfMan->pPars->fDelayOptLut && !pIfMan->pPars->fDsdBalance && !pIfMan->pPars->fUseTtPerm &&
         !pIfMan->pPars->pLutStruct && !pIfMan->pPars->fUserRecLib && !pIfMan->pPars->fUserSesLib && !pIfMan->pPars->fUserLutDec && !pIfMan->pPars->nGateSize )
        If_CutRotatePins( pIfMan, pCutBest );


    if ( pIfMan->pPars->fUseCnfs || pIfMan->pPars->fUseMv )
    {
        If_CutForEachLeafReverse( pIfMan, pCutBest, pIfLeaf, i )
            Abc_ObjAddFanin( pNodeNew, Abc_NodeFromIfDec_rec(pNtk, pNtkNew, pIfMan, pIfLeaf, vCover) );
    }
    else
    {
        If_CutForEachLeaf( pIfMan, pCutBest, pIfLeaf, i )
            Abc_ObjAddFanin( pNodeNew, Abc_NodeFromIfDec_rec(pNtk, pNtkNew, pIfMan, pIfLeaf, vCover) );
    }
    // set the level of the new node
    pNodeNew->Level = Abc_ObjLevelNew( pNodeNew );
    // derive the function of this node
    if ( pIfMan->pPars->fTruth )
    {
        if ( pIfMan->pPars->fUseBdds )
        {
            // transform truth table into the BDD
#ifdef ABC_USE_CUDD
            pNodeNew->pData = Kit_TruthToBdd( (DdManager *)pNtkNew->pManFunc, If_CutTruth(pIfMan, pCutBest), If_CutLeaveNum(pCutBest), 0 );  Cudd_Ref((DdNode *)pNodeNew->pData);
#endif
        }
        else if ( pIfMan->pPars->fUseCnfs || pIfMan->pPars->fUseMv )
        {
            // transform truth table into the BDD
#ifdef ABC_USE_CUDD
            pNodeNew->pData = Kit_TruthToBdd( (DdManager *)pNtkNew->pManFunc, If_CutTruth(pIfMan, pCutBest), If_CutLeaveNum(pCutBest), 1 );  Cudd_Ref((DdNode *)pNodeNew->pData);
#endif
        }
        else if ( pIfMan->pPars->fUseSops || pIfMan->pPars->nGateSize > 0 )
        {
            // transform truth table into the SOP
            int RetValue = Kit_TruthIsop( If_CutTruth(pIfMan, pCutBest), If_CutLeaveNum(pCutBest), vCover, 1 );
            assert( RetValue == 0 || RetValue == 1 );
            // check the case of constant cover
            if ( Vec_IntSize(vCover) == 0 || (Vec_IntSize(vCover) == 1 && Vec_IntEntry(vCover,0) == 0) )
            {
                assert( RetValue == 0 );
                pNodeNew->pData = Abc_SopCreateAnd( (Mem_Flex_t *)pNtkNew->pManFunc, If_CutLeaveNum(pCutBest), NULL );
                pNodeNew = (Vec_IntSize(vCover) == 0) ? Abc_NtkCreateNodeConst0(pNtkNew) : Abc_NtkCreateNodeConst1(pNtkNew);
            }
            else
            {
                // derive the AIG for that tree
                pNodeNew->pData = Abc_SopCreateFromIsop( (Mem_Flex_t *)pNtkNew->pManFunc, If_CutLeaveNum(pCutBest), vCover );
                if ( RetValue )
                    Abc_SopComplement( (char *)pNodeNew->pData );
            }
        }
        else if ( pIfMan->pPars->fDelayOpt )
            pNodeNew->pData = Abc_NodeBuildFromMini( (Hop_Man_t *)pNtkNew->pManFunc, pIfMan, pCutBest, 0 );
        else if ( pIfMan->pPars->fDsdBalance )
            pNodeNew->pData = Abc_NodeBuildFromMini( (Hop_Man_t *)pNtkNew->pManFunc, pIfMan, pCutBest, 1 );
        else if ( pIfMan->pPars->fUserRecLib )
        {
            extern Hop_Obj_t * Abc_RecToHop3( Hop_Man_t * pMan, If_Man_t * pIfMan, If_Cut_t * pCut, If_Obj_t * pIfObj );
            pNodeNew->pData = Abc_RecToHop3( (Hop_Man_t *)pNtkNew->pManFunc, pIfMan, pCutBest, pIfObj );
        }
        else if ( pIfMan->pPars->fUserLutDec )
        {
            extern Hop_Obj_t * Abc_DecRecordToHopCros( Abc_Ntk_t * pNtk, Hop_Man_t * pMan, If_Man_t * pIfMan, If_Cut_t * pCut, If_Obj_t * pIfObj, Vec_Int_t * vMemory );
            pNodeNew->pData = Abc_DecRecordToHopCros( pNtk, (Hop_Man_t *)pNtkNew->pManFunc, pIfMan, pCutBest, pIfObj, vCover );
        }
        else
        {
            extern Hop_Obj_t * Kit_TruthToHop( Hop_Man_t * pMan, unsigned * pTruth, int nVars, Vec_Int_t * vMemory );
            word * pTruth = If_CutTruthW(pIfMan, pCutBest);
            if ( pIfMan->pPars->fUseTtPerm )
                for ( i = 0; i < (int)pCutBest->nLeaves; i++ )
                    if ( If_CutLeafBit(pCutBest, i) )
                        Abc_TtFlip( pTruth, Abc_TtWordNum(pCutBest->nLeaves), i );
            pNodeNew->pData = Kit_TruthToHop( (Hop_Man_t *)pNtkNew->pManFunc, (unsigned *)pTruth, If_CutLeaveNum(pCutBest), vCover );
//            if ( pIfMan->pPars->fUseBat )
//                Bat_ManFuncPrintCell( *pTruth );
        }
        // complement the node if the cut was complemented
        if ( pCutBest->fCompl && !pIfMan->pPars->fDelayOpt && !pIfMan->pPars->fDsdBalance )
            Abc_NodeComplement( pNodeNew );
    }
    else
    {
        pNodeNew->pData = Abc_NodeIfToHop( (Hop_Man_t *)pNtkNew->pManFunc, pIfMan, pIfObj );
    }
    If_ObjSetCopy( pIfObj, pNodeNew );
    return pNodeNew;
}

/**Function*************************************************************

  Synopsis    [Derive one node after FPGA mapping.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
Abc_Obj_t * Abc_NodeFromIf_rec( Abc_Ntk_t * pNtkNew, If_Man_t * pIfMan, If_Obj_t * pIfObj, Vec_Int_t * vCover )
{
    Abc_Obj_t * pNodeNew;
    If_Cut_t * pCutBest;
    If_Obj_t * pIfLeaf;
    int i;
    // return if the result if known
    pNodeNew = (Abc_Obj_t *)If_ObjCopy( pIfObj );
    if ( pNodeNew )
        return pNodeNew;
    assert( pIfObj->Type == IF_AND );
    /*printf("ObjId: %i\n", pIfObj->Id);
    // Create a subnetwork for an pObj
    If_DecSubNtk_t *subnetwork = ABC_ALLOC(If_DecSubNtk_t, 1);
    subnetwork->id = pIfObj->Id;
    subnetwork->luts = Vec_PtrAlloc(10); // Or whatever size you need.
    // You can push the LUTs to the subnetwork like this:
    If_DecLut_t *lut = ABC_ALLOC(If_DecLut_t, 1);
    lut->truth_table = 1;
    lut->fan_in_luts = Vec_PtrAlloc(10);
    Vec_PtrPush(subnetwork->luts, lut);
    // Free every lut in every subnetwork
    Vec_PtrFree( lut->fan_in_luts );
    ABC_FREE( lut );
    Vec_PtrFree( subnetwork->luts );
    ABC_FREE( subnetwork );*/

    // get the parameters of the best cut
    pCutBest = If_ObjCutBest( pIfObj );
    if ( pIfMan->pPars->fUserSesLib )
    {
        // create the subgraph composed of Abc_Obj_t nodes based on the given cut
        Abc_Obj_t * pFanins[IF_MAX_FUNC_LUTSIZE];
        If_CutForEachLeaf( pIfMan, pCutBest, pIfLeaf, i )
            pFanins[i] = Abc_NodeFromIf_rec(pNtkNew, pIfMan, pIfLeaf, vCover);
        pNodeNew = Abc_ExactBuildNode( If_CutTruthW(pIfMan, pCutBest), If_CutLeaveNum(pCutBest), If_CutArrTimeProfile(pIfMan, pCutBest), pFanins, pNtkNew );
        If_ObjSetCopy( pIfObj, pNodeNew );
        return pNodeNew;
    }
    //
    // If_DecSubNtk_t* subnetwork = ABC_ALLOC(If_DecSubNtk_t, 1);
    // pIfObj->vVirtualLut = subnetwork;
    // pIfObj->vVirtualLut = NULL;
    // ABC_FREE(subnetwork);
    // create a new node
    // Instead of creating one new node: Crate a subnetwork by creating nodes and their fan-ins
    // The number of created nodes is the size of vVirtual + 1 ( for the root node )
    // the 'virtual'(VirtualLUT) and 'real'(IfObj) Fanins are stored inside the new struct VirtualLUT
    pNodeNew = Abc_NtkCreateNode( pNtkNew );
//    if ( pIfMan->pPars->pLutLib && pIfMan->pPars->pLutLib->fVarPinDelays )
    if ( !pIfMan->pPars->fDelayOpt && !pIfMan->pPars->fDelayOptLut && !pIfMan->pPars->fDsdBalance && !pIfMan->pPars->fUseTtPerm && 
         !pIfMan->pPars->pLutStruct && !pIfMan->pPars->fUserRecLib && !pIfMan->pPars->fUserSesLib && !pIfMan->pPars->fUserLutDec && !pIfMan->pPars->nGateSize )
        If_CutRotatePins( pIfMan, pCutBest );
    if ( pIfMan->pPars->fUseCnfs || pIfMan->pPars->fUseMv )
    {
        If_CutForEachLeafReverse( pIfMan, pCutBest, pIfLeaf, i )
            Abc_ObjAddFanin( pNodeNew, Abc_NodeFromIf_rec(pNtkNew, pIfMan, pIfLeaf, vCover) );
    }
    else
    {
        // ForEachLeaf of the original cut the right pNodeNew has to be selected to get his fan-in created
        If_CutForEachLeaf( pIfMan, pCutBest, pIfLeaf, i )
            Abc_ObjAddFanin( pNodeNew, Abc_NodeFromIf_rec(pNtkNew, pIfMan, pIfLeaf, vCover) );
    }
    // set the level of the new node
    pNodeNew->Level = Abc_ObjLevelNew( pNodeNew );
    // derive the function of this node
    if ( pIfMan->pPars->fTruth )
    {
        if ( pIfMan->pPars->fUseBdds )
        { 
            // transform truth table into the BDD 
#ifdef ABC_USE_CUDD
            pNodeNew->pData = Kit_TruthToBdd( (DdManager *)pNtkNew->pManFunc, If_CutTruth(pIfMan, pCutBest), If_CutLeaveNum(pCutBest), 0 );  Cudd_Ref((DdNode *)pNodeNew->pData); 
#endif
        }
        else if ( pIfMan->pPars->fUseCnfs || pIfMan->pPars->fUseMv )
        { 
            // transform truth table into the BDD 
#ifdef ABC_USE_CUDD
            pNodeNew->pData = Kit_TruthToBdd( (DdManager *)pNtkNew->pManFunc, If_CutTruth(pIfMan, pCutBest), If_CutLeaveNum(pCutBest), 1 );  Cudd_Ref((DdNode *)pNodeNew->pData); 
#endif
        }
        else if ( pIfMan->pPars->fUseSops || pIfMan->pPars->nGateSize > 0 ) 
        {
            // transform truth table into the SOP
            int RetValue = Kit_TruthIsop( If_CutTruth(pIfMan, pCutBest), If_CutLeaveNum(pCutBest), vCover, 1 );
            assert( RetValue == 0 || RetValue == 1 );
            // check the case of constant cover
            if ( Vec_IntSize(vCover) == 0 || (Vec_IntSize(vCover) == 1 && Vec_IntEntry(vCover,0) == 0) )
            {
                assert( RetValue == 0 );
                pNodeNew->pData = Abc_SopCreateAnd( (Mem_Flex_t *)pNtkNew->pManFunc, If_CutLeaveNum(pCutBest), NULL );
                pNodeNew = (Vec_IntSize(vCover) == 0) ? Abc_NtkCreateNodeConst0(pNtkNew) : Abc_NtkCreateNodeConst1(pNtkNew);
            }
            else
            {
                // derive the AIG for that tree
                pNodeNew->pData = Abc_SopCreateFromIsop( (Mem_Flex_t *)pNtkNew->pManFunc, If_CutLeaveNum(pCutBest), vCover );
                if ( RetValue )
                    Abc_SopComplement( (char *)pNodeNew->pData );
            }
        }
        else if ( pIfMan->pPars->fDelayOpt )
            pNodeNew->pData = Abc_NodeBuildFromMini( (Hop_Man_t *)pNtkNew->pManFunc, pIfMan, pCutBest, 0 );
        else if ( pIfMan->pPars->fDsdBalance )
            pNodeNew->pData = Abc_NodeBuildFromMini( (Hop_Man_t *)pNtkNew->pManFunc, pIfMan, pCutBest, 1 );
        else if ( pIfMan->pPars->fUserRecLib )
        {
            extern Hop_Obj_t * Abc_RecToHop3( Hop_Man_t * pMan, If_Man_t * pIfMan, If_Cut_t * pCut, If_Obj_t * pIfObj );
            pNodeNew->pData = Abc_RecToHop3( (Hop_Man_t *)pNtkNew->pManFunc, pIfMan, pCutBest, pIfObj ); 
        }
        else if ( pIfMan->pPars->fUserLutDec )
        {
            extern Hop_Obj_t * Abc_DecRecordToHop( Hop_Man_t * pMan, If_Man_t * pIfMan, If_Cut_t * pCut, If_Obj_t * pIfObj, Vec_Int_t * vMemory );
            pNodeNew->pData = Abc_DecRecordToHop( (Hop_Man_t *)pNtkNew->pManFunc, pIfMan, pCutBest, pIfObj, vCover ); 
        }
        else
        {
            extern Hop_Obj_t * Kit_TruthToHop( Hop_Man_t * pMan, unsigned * pTruth, int nVars, Vec_Int_t * vMemory );
            word * pTruth = If_CutTruthW(pIfMan, pCutBest);
            if ( pIfMan->pPars->fUseTtPerm )
                for ( i = 0; i < (int)pCutBest->nLeaves; i++ )
                    if ( If_CutLeafBit(pCutBest, i) )
                        Abc_TtFlip( pTruth, Abc_TtWordNum(pCutBest->nLeaves), i );
            pNodeNew->pData = Kit_TruthToHop( (Hop_Man_t *)pNtkNew->pManFunc, (unsigned *)pTruth, If_CutLeaveNum(pCutBest), vCover );
//            if ( pIfMan->pPars->fUseBat )
//                Bat_ManFuncPrintCell( *pTruth );
        }
        // complement the node if the cut was complemented
        if ( pCutBest->fCompl && !pIfMan->pPars->fDelayOpt && !pIfMan->pPars->fDsdBalance )
            Abc_NodeComplement( pNodeNew );
    }
    else
    {
        pNodeNew->pData = Abc_NodeIfToHop( (Hop_Man_t *)pNtkNew->pManFunc, pIfMan, pIfObj );
    }
    If_ObjSetCopy( pIfObj, pNodeNew );
    return pNodeNew;
}

/**Function*************************************************************

  Synopsis    [Recursively derives the truth table for the cut.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
Hop_Obj_t * Abc_NodeIfToHop_rec( Hop_Man_t * pHopMan, If_Man_t * pIfMan, If_Obj_t * pIfObj, Vec_Ptr_t * vVisited )
{
    If_Cut_t * pCut;
    Hop_Obj_t * gFunc, * gFunc0, * gFunc1;
    // get the best cut
    pCut = If_ObjCutBest(pIfObj);
    // if the cut is visited, return the result
    if ( If_CutData(pCut) )
        return (Hop_Obj_t *)If_CutData(pCut);
    // compute the functions of the children
    gFunc0 = Abc_NodeIfToHop_rec( pHopMan, pIfMan, pIfObj->pFanin0, vVisited );
    gFunc1 = Abc_NodeIfToHop_rec( pHopMan, pIfMan, pIfObj->pFanin1, vVisited );
    // get the function of the cut
    gFunc  = Hop_And( pHopMan, Hop_NotCond(gFunc0, pIfObj->fCompl0), Hop_NotCond(gFunc1, pIfObj->fCompl1) );  
    assert( If_CutData(pCut) == NULL );
    If_CutSetData( pCut, gFunc );
    // add this cut to the visited list
    Vec_PtrPush( vVisited, pCut );
    return gFunc;
}


/**Function*************************************************************

  Synopsis    [Recursively derives the truth table for the cut.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
Hop_Obj_t * Abc_NodeIfToHop2_rec( Hop_Man_t * pHopMan, If_Man_t * pIfMan, If_Obj_t * pIfObj, Vec_Ptr_t * vVisited )
{
    If_Cut_t * pCut;
    If_Obj_t * pTemp;
    Hop_Obj_t * gFunc, * gFunc0, * gFunc1;
    // get the best cut
    pCut = If_ObjCutBest(pIfObj);
    // if the cut is visited, return the result
    if ( If_CutData(pCut) )
        return (Hop_Obj_t *)If_CutData(pCut);
    // mark the node as visited
    Vec_PtrPush( vVisited, pCut );
    // insert the worst case
    If_CutSetData( pCut, (void *)1 );
    // skip in case of primary input
    if ( If_ObjIsCi(pIfObj) )
        return (Hop_Obj_t *)If_CutData(pCut);
    // compute the functions of the children
    for ( pTemp = pIfObj; pTemp; pTemp = pTemp->pEquiv )
    {
        gFunc0 = Abc_NodeIfToHop2_rec( pHopMan, pIfMan, pTemp->pFanin0, vVisited );
        if ( gFunc0 == (void *)1 )
            continue;
        gFunc1 = Abc_NodeIfToHop2_rec( pHopMan, pIfMan, pTemp->pFanin1, vVisited );
        if ( gFunc1 == (void *)1 )
            continue;
        // both branches are solved
        gFunc = Hop_And( pHopMan, Hop_NotCond(gFunc0, pTemp->fCompl0), Hop_NotCond(gFunc1, pTemp->fCompl1) ); 
        if ( pTemp->fPhase != pIfObj->fPhase )
            gFunc = Hop_Not(gFunc);
        If_CutSetData( pCut, gFunc );
        break;
    }
    return (Hop_Obj_t *)If_CutData(pCut);
}

/**Function*************************************************************

  Synopsis    [Derives the truth table for one cut.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
Hop_Obj_t * Abc_NodeIfToHop( Hop_Man_t * pHopMan, If_Man_t * pIfMan, If_Obj_t * pIfObj )
{
    If_Cut_t * pCut;
    Hop_Obj_t * gFunc;
    If_Obj_t * pLeaf;
    int i;
    // get the best cut
    pCut = If_ObjCutBest(pIfObj);
    assert( pCut->nLeaves > 1 );
    // set the leaf variables
    If_CutForEachLeaf( pIfMan, pCut, pLeaf, i )
        If_CutSetData( If_ObjCutBest(pLeaf), Hop_IthVar(pHopMan, i) );
    // recursively compute the function while collecting visited cuts
    Vec_PtrClear( pIfMan->vTemp );
    gFunc = Abc_NodeIfToHop2_rec( pHopMan, pIfMan, pIfObj, pIfMan->vTemp ); 
    if ( gFunc == (void *)1 )
    {
        printf( "Abc_NodeIfToHop(): Computing local AIG has failed.\n" );
        return NULL;
    }
    // clean the cuts
    If_CutForEachLeaf( pIfMan, pCut, pLeaf, i )
        If_CutSetData( If_ObjCutBest(pLeaf), NULL );
    Vec_PtrForEachEntry( If_Cut_t *, pIfMan->vTemp, pCut, i )
        If_CutSetData( pCut, NULL );
    return gFunc;
}


/**Function*************************************************************

  Synopsis    [Comparison for two nodes with the flow.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int Abc_ObjCompareFlow( Abc_Obj_t ** ppNode0, Abc_Obj_t ** ppNode1 )
{
    float Flow0 = Abc_Int2Float((int)(ABC_PTRINT_T)(*ppNode0)->pCopy);
    float Flow1 = Abc_Int2Float((int)(ABC_PTRINT_T)(*ppNode1)->pCopy);
    if ( Flow0 > Flow1 )
        return -1;
    if ( Flow0 < Flow1 )
        return 1;
    return 0;
}

/**Function*************************************************************

  Synopsis    [Orders AIG nodes so that nodes from larger cones go first.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Abc_NtkFindGoodOrder_rec( Abc_Obj_t * pNode, Vec_Ptr_t * vNodes )
{
    if ( !Abc_ObjIsNode(pNode) )
        return;
    assert( Abc_ObjIsNode( pNode ) );
    // if this node is already visited, skip
    if ( Abc_NodeIsTravIdCurrent( pNode ) )
        return;
    // mark the node as visited
    Abc_NodeSetTravIdCurrent( pNode );
    // visit the transitive fanin of the node
    Abc_NtkFindGoodOrder_rec( Abc_ObjFanin0(pNode), vNodes );
    Abc_NtkFindGoodOrder_rec( Abc_ObjFanin1(pNode), vNodes );
    // add the node after the fanins have been added
    Vec_PtrPush( vNodes, pNode );
}

/**Function*************************************************************

  Synopsis    [Orders AIG nodes so that nodes from larger cones go first.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
Vec_Ptr_t * Abc_NtkFindGoodOrder( Abc_Ntk_t * pNtk )
{
    Vec_Ptr_t * vNodes, * vCos;
    Abc_Obj_t * pNode, * pFanin0, * pFanin1;
    float Flow0, Flow1;
    int i;

    // initialize the flow
    Abc_AigConst1(pNtk)->pCopy = NULL;
    Abc_NtkForEachCi( pNtk, pNode, i )
        pNode->pCopy = NULL;
    // compute the flow
    Abc_AigForEachAnd( pNtk, pNode, i )
    {
        pFanin0 = Abc_ObjFanin0(pNode);
        pFanin1 = Abc_ObjFanin1(pNode);
        Flow0 = Abc_Int2Float((int)(ABC_PTRINT_T)pFanin0->pCopy)/Abc_ObjFanoutNum(pFanin0);
        Flow1 = Abc_Int2Float((int)(ABC_PTRINT_T)pFanin1->pCopy)/Abc_ObjFanoutNum(pFanin1);
        pNode->pCopy = (Abc_Obj_t *)(ABC_PTRINT_T)Abc_Float2Int(Flow0 + Flow1+(float)1.0);
    }
    // find the flow of the COs
    vCos = Vec_PtrAlloc( Abc_NtkCoNum(pNtk) );
    Abc_NtkForEachCo( pNtk, pNode, i )
    {
        pNode->pCopy = Abc_ObjFanin0(pNode)->pCopy;
//        pNode->pCopy = (Abc_Obj_t *)Abc_Float2Int((float)Abc_ObjFanin0(pNode)->Level);
        Vec_PtrPush( vCos, pNode );
    }

    // sort nodes in the increasing order of the flow
    qsort( (Abc_Obj_t **)Vec_PtrArray(vCos), (size_t)Abc_NtkCoNum(pNtk), 
        sizeof(Abc_Obj_t *), (int (*)(const void *, const void *))Abc_ObjCompareFlow );
    // verify sorting
    pFanin0 = (Abc_Obj_t *)Vec_PtrEntry(vCos, 0);
    pFanin1 = (Abc_Obj_t *)Vec_PtrEntryLast(vCos);
    assert( Abc_Int2Float((int)(ABC_PTRINT_T)pFanin0->pCopy) >= Abc_Int2Float((int)(ABC_PTRINT_T)pFanin1->pCopy) );

    // collect the nodes in the topological order from the new array
    Abc_NtkIncrementTravId( pNtk );
    vNodes = Vec_PtrAlloc( 100 );
    Vec_PtrForEachEntry( Abc_Obj_t *, vCos, pNode, i )
        Abc_NtkFindGoodOrder_rec( Abc_ObjFanin0(pNode), vNodes );
    Vec_PtrFree( vCos );
    return vNodes;
}


/**Function*************************************************************

  Synopsis    [Sets PO drivers.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Abc_NtkMarkMux( Abc_Obj_t * pDriver, Abc_Obj_t ** ppNode1, Abc_Obj_t ** ppNode2 )
{
    Abc_Obj_t * pNodeC, * pNodeT, * pNodeE;
    If_Obj_t * pIfObj;

    *ppNode1 = NULL;
    *ppNode2 = NULL;
    if ( pDriver == NULL )
        return;
    if ( !Abc_NodeIsMuxType(pDriver) )
        return;

    pNodeC = Abc_NodeRecognizeMux( pDriver, &pNodeT, &pNodeE );

    pIfObj = If_Regular( (If_Obj_t *)Abc_ObjFanin0(pDriver)->pCopy );
    if ( If_ObjIsAnd(pIfObj) )
        pIfObj->fSkipCut = 1;
    pIfObj = If_Regular( (If_Obj_t *)Abc_ObjFanin1(pDriver)->pCopy );
    if ( If_ObjIsAnd(pIfObj) )
        pIfObj->fSkipCut = 1;

    pIfObj = If_Regular( (If_Obj_t *)Abc_ObjRegular(pNodeC)->pCopy );
    if ( If_ObjIsAnd(pIfObj) )
        pIfObj->fSkipCut = 1;

/*
    pIfObj = If_Regular( (If_Obj_t *)Abc_ObjRegular(pNodeT)->pCopy );
    if ( If_ObjIsAnd(pIfObj) )
        pIfObj->fSkipCut = 1;
    pIfObj = If_Regular( (If_Obj_t *)Abc_ObjRegular(pNodeE)->pCopy );
    if ( If_ObjIsAnd(pIfObj) )
        pIfObj->fSkipCut = 1;
*/
    *ppNode1 = Abc_ObjRegular(pNodeC);
    *ppNode2 = Abc_ObjRegular(pNodeT);
}


////////////////////////////////////////////////////////////////////////
///                       END OF FILE                                ///
////////////////////////////////////////////////////////////////////////


ABC_NAMESPACE_IMPL_END

