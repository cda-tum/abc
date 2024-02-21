//
// Created by benjamin on 2/12/24.
//

#ifndef ABC_IFCROSDEC_H
#define ABC_IFCROSDEC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

extern If_Man_t *  Abc_NtkToIf( Abc_Ntk_t * pNtk, If_Par_t * pPars );
extern Abc_Obj_t * Abc_NodeFromIf_rec( Abc_Ntk_t * pNtkNew, If_Man_t * pIfMan, If_Obj_t * pIfObj, Vec_Int_t * vCover );
static Hop_Obj_t * Abc_NodeIfToHop( Hop_Man_t * pHopMan, If_Man_t * pIfMan, If_Obj_t * pIfObj );
extern If_DecSubNtk_t * Abc_DecRecordToHopCros4( Abc_Ntk_t * pNtk, Hop_Man_t * pMan, If_Man_t * pIfMan, If_Cut_t * pCutBest, If_Obj_t * pIfObj, If_Obj_t * pIfObjCur, Vec_Int_t * vCover, Vec_Int_t * vUsedNodes );
extern If_DecSubNtk_t * Abc_DecRecordToHopCros3( Abc_Ntk_t * pNtk, Hop_Man_t * pMan, If_Man_t * pIfMan, If_Cut_t * pCutBest, If_Obj_t * pIfObj, If_Obj_t * pIfObjCur, Vec_Int_t * vCover, Vec_Int_t * vUsedNodes );

extern void Abc_NtkBddReorder( Abc_Ntk_t * pNtk, int fVerbose );
extern void Abc_NtkBidecResyn( Abc_Ntk_t * pNtk, int fVerbose );

/*typedef struct {
    word* sequence;
    int size;
} SequenceData;
static SequenceData global3InpSequence;
void generateAndStoreBinarySequence(int num_bits) {
    global3InpSequence.sequence = generateBinarySequence(num_bits, &global3InpSequence.size);
}
generateAndStoreBinarySequence(3);*/
/**Function*************************************************************

  Synopsis    []

  Description [Used helper functions for transition.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
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
If_DecSubNtk_t * Abc_TraverseNodesIf_rec( Abc_Ntk_t * pNtk, If_Man_t * pIfMan, If_Obj_t * pIfObj, If_Obj_t * pIfObjStart, If_Cut_t * pC0,
                                          Hop_Man_t * pMan, Vec_Int_t * vCover, Vec_Int_t * vUsedNodes )
{
    If_Cut_t * pC1;
    If_Obj_t * pIfLeaf;
    If_DecSubNtk_t * pDec;
    int i;

    if( (pIfObj->Type != IF_AND) || pIfObj->fVisit == 1 )
        return NULL;

    pC1 = If_ObjCutBest( pIfObj );

    If_CutForEachLeaf(pIfMan, pC1, pIfLeaf, i)
    {
        pDec = Abc_TraverseNodesIf_rec( pNtk, pIfMan, pIfLeaf, pIfObjStart, pC0, pMan, vCover, vUsedNodes ); // recursive call
        if (pDec != NULL && pC0 != pC1) {
            return pDec; // break and return if a match is found
        }
    }

    int sharedLeaves = CountSharedLeaves( pC0, pC1 );
    if( sharedLeaves == 4 && pC1 != pC0 )
    {
        assert( pIfObj->Type == IF_AND );
        pDec = Abc_DecRecordToHopCros4( pNtk, pMan, pIfMan, pC0, pIfObjStart, pIfObj, vCover, vUsedNodes );
        // printf("Number of shared leaves: %d\n", sharedLeaves);
        return pDec;
    }
    /*if( sharedLeaves == 3 && pC1 != pC0 )
    {
        assert( pIfObj->Type == IF_AND );
        pDec = Abc_DecRecordToHopCros3( pNtk, pMan, pIfMan, pC0, pIfObjStart, pIfObj, vCover, vUsedNodes );
        // printf("Number of shared leaves: %d\n", sharedLeaves);
        return pDec;
    }*/

    pIfObj->fVisit = 1;
    // Last cut can be skipped since there is no cut to compare it to
    if ( pC1 == If_ObjCutBest( If_ObjFanin0( If_ManCo( pIfMan, Vec_PtrSize( pIfMan->vCos) - 1 ) ) ) )
    {
        return NULL;
    }

    return NULL;
}
/**Function*************************************************************

  Synopsis    []

  Description [Used helper functions for bit operations.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void printBits(word input) {
    int bits_in_word = 8 * sizeof(word);
    for (int bit = 0; bit < bits_in_word; ++bit) {
        printf("%lu", (input >> (bits_in_word - bit - 1)) & 1);
    }
    printf("\n");
}
word convertTo64Bit(word smallTruth, int bits) {
    word result = 0;
    int repeats = 64 / bits;

    for (int i = 0; i < repeats; i++)
        result = (result << bits) | smallTruth;

    return result;
}
// Function to generate the binary sequence and return it as an array
word* generateBinarySequence(int num_bits, int* size) {
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
            int indices[] = {2, 1};
            swapAdjacentIterations(pTruth, nWords, indices, sizeof(indices) / sizeof(indices[0]));
            break;
        }
        case 3: {
            int indices[] = {2, 1, 0};
            swapAdjacentIterations(pTruth, nWords, indices, sizeof(indices) / sizeof(indices[0]));
            break;
        }

        default:
            printf("Invalid bit value %i\n", nShifts);
            break;
    }
}
void calculate_bit_nums_2_(int iter, int* bit_num0, int* bit_num1) {
    switch(iter) {
        case 0:
            *bit_num0 = 2;
            *bit_num1 = 3;
            break;
        case 1:
            *bit_num0 = 1;
            *bit_num1 = 3;
            break;
        case 2:
            *bit_num0 = 1;
            *bit_num1 = 2;
            break;
        case 3:
            *bit_num0 = 0;
            *bit_num1 = 3;
            break;
        case 4:
            *bit_num0 = 0;
            *bit_num1 = 2;
            break;
        case 5:
            *bit_num0 = 0;
            *bit_num1 = 1;
            break;
        default:
            printf("Invalid iter value. Must be between 0 and 5.\n");
            break;
    }
}
int apply_operations_3_2( word lut4inp, word lut3inp, word lut2inp, int iter )
{
    int bit_num = 3 - iter;
    word bit_mask = s_Truths6[bit_num];
    /*word bit1 = s_Truths6[3];
    word bit2 = s_Truths6[2] ;
    word bit3 = s_Truths6[1] ;
    word bit4 = s_Truths6[0] ;*/

    // bit1 = convertTo64Bit(bit1, 16);
    word output = 0;
    int bits_in_word = 4;
    for (int bit = 0; bit < bits_in_word; ++bit)
    {
        int do_operation = (int)(lut2inp >> (bits_in_word - bit - 1)) & 1;
        if ( do_operation )
        {
            switch (bit)
            {
                case 0:
                    output = output | ( bit_mask & lut3inp );
                    break;

                case 1:
                    output = output | ( bit_mask & ~lut3inp );
                    break;

                case 2:
                    output = output | ( ~bit_mask & lut3inp );
                    break;

                case 3:
                    output = output | ( ~bit_mask & ~lut3inp );
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
int apply_operations_3_3( word lut4inp, word lut3inp0, word lut3inp1, int iter )
{
    int bit_num = 3 - iter;
    word bit_mask0, bit_mask1;
    bit_mask0 = s_Truths6[bit_num];

    word output;
    for (int i = 0; i < 4; ++i)
    {
        output = 0;
        if ( i == bit_num )
            continue;
        bit_mask1 = s_Truths6[i];
        int bits_in_word = 8;  // Bits in word now 8 to accommodate 3-input LUT
        for (int bit = 0; bit < bits_in_word; ++bit)
        {
            int do_operation = (int)(lut3inp1 >> (bits_in_word - bit - 1)) & 1;
            if ( do_operation )
            {
                switch (bit)
                {
                    case 0:
                        output = output | ( bit_mask0 & bit_mask1 & lut3inp0 );
                        break;
                    case 1:
                        output = output | ( bit_mask0 & bit_mask1 & ~lut3inp0 );
                        break;
                    case 2:
                        output = output | ( bit_mask0 & ~bit_mask1 & lut3inp0 );
                        break;
                    case 3:
                        output = output | ( bit_mask0 & ~bit_mask1 & ~lut3inp0 );
                        break;
                    case 4:
                        output = output | ( ~bit_mask0 & bit_mask1 & lut3inp0 );
                        break;
                    case 5:
                        output = output | ( ~bit_mask0 & bit_mask1 & ~lut3inp0 );
                        break;
                    case 6:
                        output = output | ( ~bit_mask0 & ~bit_mask1 & lut3inp0 );
                        break;
                    case 7:
                        output = output | ( ~bit_mask0 & ~bit_mask1 & ~lut3inp0 );
                        break;
                    default:
                        printf("Invalid bit value\n");
                        break;
                }
            }
        }
        if( lut4inp == output )
        {
            return i+1;
        }
    }
    return 0;
}
int apply_operations_2_3( word lut4inp, word lut3inp0, word lut3inp1, int iter )
{
    // the iteration here determines both other bit masks
    int bit_num0, bit_num1;
    calculate_bit_nums_2_( iter, &bit_num0, &bit_num1 );
    word bit_mask0, bit_mask1;
    bit_mask0 = s_Truths6[bit_num0];
    bit_mask1 = s_Truths6[bit_num1];

    word output = 0;
    int bits_in_word = 8;  // Bits in word now 8 to accommodate 3-input LUT
    for (int bit = 0; bit < bits_in_word; ++bit)
    {
        int do_operation = (int)(lut3inp1 >> (bits_in_word - bit - 1)) & 1;
        if ( do_operation )
        {
            switch (bit)
            {
                case 0:
                    output = output | ( bit_mask0 & bit_mask1 & lut3inp0 );
                    break;
                case 1:
                    output = output | ( bit_mask0 & bit_mask1 & ~lut3inp0 );
                    break;
                case 2:
                    output = output | ( bit_mask0 & ~bit_mask1 & lut3inp0 );
                    break;
                case 3:
                    output = output | ( bit_mask0 & ~bit_mask1 & ~lut3inp0 );
                    break;
                case 4:
                    output = output | ( ~bit_mask0 & bit_mask1 & lut3inp0 );
                    break;
                case 5:
                    output = output | ( ~bit_mask0 & bit_mask1 & ~lut3inp0 );
                    break;
                case 6:
                    output = output | ( ~bit_mask0 & ~bit_mask1 & lut3inp0 );
                    break;
                case 7:
                    output = output | ( ~bit_mask0 & ~bit_mask1 & ~lut3inp0 );
                    break;
                default:
                    printf("Invalid bit value\n");
                    break;
            }
        }
    if( lut4inp == output )
    {
        return 1;
    }
    }
    return 0;
}
int apply_operations_3_4( word lut4inp, word lut3inp0, word lut3inp1, int iter )
{
    int bit_num = 3 - iter;
    word bit_mask0, bit_mask1, bit_mask2;
    bit_mask0 = s_Truths6[bit_num];
    word output;

    for (int i = 0; i < 4; ++i)
    {
        if ( i == bit_num )
            continue;
        bit_mask1 = s_Truths6[i];
        for (int j = 0; j < 4; ++j)
        {
            output = 0;
            if ( j == bit_num || i == j )
                continue;
            bit_mask2 = s_Truths6[j];
            int bits_in_word = 16;  // Bits in word now 8 to accommodate 3-input LUT
            for (int bit = 0; bit < bits_in_word; ++bit)
            {
                int do_operation = (int) (lut3inp1 >> (bits_in_word - bit - 1)) & 1;
                if (do_operation)
                {
                    switch (bit)
                    {
                        case 0:
                            output = output | (bit_mask0 & bit_mask1 & bit_mask2 & lut3inp0);
                            break;
                        case 1:
                            output = output | (bit_mask0 & bit_mask1 & bit_mask2 & ~lut3inp0);
                            break;
                        case 2:
                            output = output | (bit_mask0 & bit_mask1 & ~bit_mask2 & lut3inp0);
                            break;
                        case 3:
                            output = output | (bit_mask0 & bit_mask1 & ~bit_mask2 & ~lut3inp0);
                            break;
                        case 4:
                            output = output | (bit_mask0 & ~bit_mask1 & bit_mask2 & lut3inp0);
                            break;
                        case 5:
                            output = output | (bit_mask0 & ~bit_mask1 & bit_mask2 & ~lut3inp0);
                            break;
                        case 6:
                            output = output | (bit_mask0 & ~bit_mask1 & ~bit_mask2 & lut3inp0);
                            break;
                        case 7:
                            output = output | (bit_mask0 & ~bit_mask1 & ~bit_mask2 & ~lut3inp0);
                            break;
                        case 8:
                            output = output | (~bit_mask0 & bit_mask1 & bit_mask2 & lut3inp0);
                            break;
                        case 9:
                            output = output | (~bit_mask0 & bit_mask1 & bit_mask2 & ~lut3inp0);
                            break;
                        case 10:
                            output = output | (~bit_mask0 & bit_mask1 & ~bit_mask2 & lut3inp0);
                            break;
                        case 11:
                            output = output | (~bit_mask0 & bit_mask1 & ~bit_mask2 & ~lut3inp0);
                            break;
                        case 12:
                            output = output | (~bit_mask0 & ~bit_mask1 & bit_mask2 & lut3inp0);
                            break;
                        case 13:
                            output = output | (~bit_mask0 & ~bit_mask1 & bit_mask2 & ~lut3inp0);
                            break;
                        case 14:
                            output = output | (~bit_mask0 & ~bit_mask1 & ~bit_mask2 & lut3inp0);
                            break;
                        case 15:
                            output = output | (~bit_mask0 & ~bit_mask1 & ~bit_mask2 & ~lut3inp0);
                            break;
                        default:
                            printf("Invalid case value\n");
                            break;
                    }
                }
            }
            if (lut4inp == output)
            {
                return i*10 + j + 1;
            }
        }
    }
    return 0;
}
/**Function*************************************************************

  Synopsis    []

  Description [Used helper functions for memory management.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void freeDecNtk(If_DecSubNtk_t* DecNtk) {
    if(DecNtk)
    {
        if(DecNtk->vLut)
        {
            Vec_PtrFree(DecNtk->vLut);
        }
        free(DecNtk);
    }
}
void freeSubnetwork(If_DecSubNtk_t* subNtk) {
    if (subNtk->vLut)
    {
        for (int k = 0; k < Vec_PtrSize(subNtk->vLut); k++)
        {
            If_DecLut_t* lut = (If_DecLut_t*) Vec_PtrEntry(subNtk->vLut, k);
            if(lut)
            {
                free(lut->vRealFis);
                free(lut);
            }
        }
        Vec_PtrFree(subNtk->vLut);
    }
    free(subNtk);
}
void freeSubnetworks(Vec_Ptr_t* vDecSubNtk) {
    printf( "%i subnetworks freed\n", Vec_PtrSize(vDecSubNtk) );

    for (int j = 0; j < Vec_PtrSize(vDecSubNtk); j++)
    {
        If_DecSubNtk_t* subNtk = (If_DecSubNtk_t*) Vec_PtrEntry(vDecSubNtk, j);
        // freeing vLut's entries, assuming it holds pointers to If_DecLut_t
        if (subNtk->vLut)
        {
            for (int k = 0; k < Vec_PtrSize(subNtk->vLut); k++)
            {
                If_DecLut_t* lut = (If_DecLut_t*) Vec_PtrEntry(subNtk->vLut, k);
                if(lut)
                {
                    free(lut->vRealFis);
                    free(lut);
                }
            }
            Vec_PtrFree(subNtk->vLut);
        }
        free(subNtk);
    }

    Vec_PtrFree(vDecSubNtk);
}
If_DecLut_t* createAndInitDecLut(If_DecSubNtk_t * DecNtk, word truth_table_value, int num_real_fis) {
    // Create a new LUT
    If_DecLut_t *pDecLut = (If_DecLut_t*)malloc(sizeof(If_DecLut_t));

    // Initialize memory
    memset(pDecLut, 0, sizeof(If_DecLut_t));

    // Assign appropriate variables
    pDecLut->truth_table = truth_table_value;

    // Allocate memory for vVirtFis and vRealFis
    pDecLut->vRealFis = malloc(num_real_fis * sizeof(int));

    Vec_PtrPush(DecNtk->vLut, pDecLut);

    return pDecLut;
}
/**Function*************************************************************

  Synopsis    []

  Description [Unused helper functions.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void If_ManCleanMarkS( If_Man_t * p )
{
    If_Obj_t * pObj;
    int i;
    If_ManForEachObj( p, pObj, i )
        pObj->fSpec = 0;
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
}
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
}*/
void printBits2(word input) {
    int bits_in_word = 4;
    for (int bit = 0; bit < bits_in_word; ++bit) {
        printf("%d", (input >> (bits_in_word - bit - 1)) & 1);
    }
    printf("\n");
}
// Function to convert 16, 8, or 4-bit sequence into 64-bit sequence by repeating it
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

/**Function*************************************************************

  Synopsis    []

  Description [Get the Hop function of the custom decomposed truth tables.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
Hop_Obj_t * Abc_DecSubNtkToHop( Hop_Man_t * pMan, If_Man_t * pIfMan, Vec_Int_t * vCover, word * pTruth, int nVars )
{
    // get the truth table
    // perform LUT-decomposition and return the LUT-structure
    // convert the LUT-structure into a set of logic nodes in Abc_Ntk_t

    // this is a placeholder, which takes the truth table and converts it into an AIG without LUT-decomposition
    extern Hop_Obj_t * Kit_TruthToHop( Hop_Man_t * pMan, unsigned * pTruth, int nVars, Vec_Int_t * vMemory );
    assert( !pIfMan->pPars->fUseTtPerm );
    return Kit_TruthToHop( (Hop_Man_t *)pMan, (unsigned *)pTruth, nVars, vCover );
}
int** reorderArrays3(const int arr[])
{
    int** reorderedArrays = (int**)malloc(4 * sizeof(int*));
    if (reorderedArrays == NULL) {
        printf("Memory allocation failed!\n");
        exit(1);
    }

    for (int i = 0; i < 4; i++) {
        reorderedArrays[i] = (int*)malloc(4 * sizeof(int));
        if (reorderedArrays[i] == NULL) {
            printf("Memory allocation failed!\n");
            exit(1);
        }
    }

    // Iterator 0: {4, 1, 2, 3}
    reorderedArrays[0][0] = arr[0]; reorderedArrays[0][1] = arr[3]; reorderedArrays[0][2] = arr[2]; reorderedArrays[0][3] = arr[1];
    // Iterator 1: {3, 1, 2, 4}
    reorderedArrays[1][0] = arr[1]; reorderedArrays[1][1] = arr[3]; reorderedArrays[1][2] = arr[2]; reorderedArrays[1][3] = arr[0];
    // Iterator 2: {2, 1, 3, 4}
    reorderedArrays[2][0] = arr[2]; reorderedArrays[2][1] = arr[3]; reorderedArrays[2][2] = arr[1]; reorderedArrays[2][3] = arr[0];
    // Iterator 3: {1, 2, 3, 4}
    reorderedArrays[3][0] = arr[3]; reorderedArrays[3][1] = arr[2]; reorderedArrays[3][2] = arr[1]; reorderedArrays[3][3] = arr[0];

    return reorderedArrays;
}

/**Function*************************************************************

  Synopsis    []

  Description [Decomposition Functions.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
static void freeSequences1(word* lut_sequence, int **ci_sequences){
    free(lut_sequence);
    for (int i = 0; i < 4; i++)
        free(ci_sequences[i]); // Free memory for each reordered array
    free(ci_sequences);
}
static void freeSequences2(word* lut_sequence0, word* lut_sequence1, int **ci_sequences){
    free(lut_sequence0);
    free(lut_sequence1);
    for (int i = 0; i < 4; i++)
        free(ci_sequences[i]); // Free memory for each reordered array
    free(ci_sequences);
}
void create_and_init_dec_lut02( If_DecSubNtk_t* DecNtk, word cur2inpLut, int is_complement, int* ci_sequence, If_DecLut_t* pDecLut2 )
{
    word tt = convertTo64Bit((word) cur2inpLut, 4);

    if (is_complement == 1 || is_complement == 3) {
        tt = ~tt;
    }

    If_DecLut_t* pDecLut0 = createAndInitDecLut(DecNtk, tt, 1);
    pDecLut0->nVars = 2;

    // Add LUT to DecSubNtk_t
    DecNtk->pLutId = pDecLut0;

    // create the right fan-in relationships
    pDecLut0->pVirtFi = pDecLut2;
    pDecLut0->vRealFis[0] = ci_sequence[0];
    pDecLut0->nRealFis = 1;
}
void create_and_init_dec_lut03(If_DecSubNtk_t* DecNtk, word cur3inpLut, int is_complement, int* ci_sequence, If_DecLut_t* pDecLut2, int* vLeaves, int found_dec)
{
    word tt = convertTo64Bit(cur3inpLut, 8);

    if (is_complement == 1 || is_complement == 3) {
        tt = ~tt;
    }

    If_DecLut_t* pDecLut0 = createAndInitDecLut(DecNtk, tt, 2);
    pDecLut0->nVars = 3;

    // Add LUT to DecSubNtk_t
    DecNtk->pLutId = pDecLut0;

    // create the right fan-in relationships
    pDecLut0->pVirtFi = pDecLut2;
    pDecLut0->vRealFis[1] = ci_sequence[0];
    pDecLut0->vRealFis[0] = vLeaves[3 - (found_dec - 1)];
    pDecLut0->nRealFis = 2;
}
void create_and_init_dec_lut12( If_DecSubNtk_t* DecNtk, word cur2inpLut, int is_complement, int* ci_sequence, If_DecLut_t* pDecLut2 )
{
    word tt = convertTo64Bit((word) cur2inpLut, 4);

    if (is_complement == 2 || is_complement == 3) {
        tt = ~tt;
    }

    If_DecLut_t* pDecLut1 = createAndInitDecLut(DecNtk, tt, 1);
    pDecLut1->nVars = 2;

    // Add LUT to DecSubNtk_t
    DecNtk->pLutSharedId = pDecLut1;

    // create the right fan-in relationships
    pDecLut1->pVirtFi = pDecLut2;
    pDecLut1->vRealFis[0] = ci_sequence[0];
    pDecLut1->nRealFis = 1;
}
void create_and_init_dec_lut13(If_DecSubNtk_t* DecNtk, word cur3inpLut, int is_complement, int* ci_sequence, If_DecLut_t* pDecLut2, int* vLeaves, int found_dec)
{
    word tt = convertTo64Bit(cur3inpLut, 8);

    if (is_complement == 2 || is_complement == 3) {
        tt = ~tt;
    }

    If_DecLut_t* pDecLut1 = createAndInitDecLut(DecNtk, tt, 2);
    pDecLut1->nVars = 3;

    // Add LUT to DecSubNtk_t
    DecNtk->pLutSharedId = pDecLut1;

    // create the right fan-in relationships
    pDecLut1->pVirtFi = pDecLut2;
    pDecLut1->vRealFis[0] = ci_sequence[0];
    pDecLut1->vRealFis[1] = vLeaves[3 - (found_dec - 1)];
    pDecLut1->nRealFis = 2;
}
void create_and_init_dec_lut14(If_DecSubNtk_t* DecNtk, word cur4inpLut, int is_complement, int* ci_sequence, If_DecLut_t* pDecLut2, int* vLeaves, int found_dec)
{
    word tt = convertTo64Bit(cur4inpLut, 16);

    if (is_complement == 2 || is_complement == 3) {
        tt = ~tt;
    }

    If_DecLut_t* pDecLut1 = createAndInitDecLut(DecNtk, tt, 3);
    pDecLut1->nVars = 4;

    // Add LUT to DecSubNtk_t
    DecNtk->pLutSharedId = pDecLut1;

    // create the right fan-in relationships
    pDecLut1->pVirtFi = pDecLut2;
    pDecLut1->vRealFis[0] = ci_sequence[0];

    // Calculate real fan-in using found_dec
    int secondDigit = found_dec % 10;
    int firstDigit = found_dec / 10;

    pDecLut1->vRealFis[1] = vLeaves[3 - firstDigit];
    pDecLut1->vRealFis[2] = vLeaves[3 - (secondDigit - 1)];
    pDecLut1->nRealFis = 3;
}
If_DecLut_t* create_and_init_dec_lut2(If_DecSubNtk_t* DecNtk, word cur3inpLut, int* ci_sequence)
{
    If_DecLut_t* pDecLut2 = createAndInitDecLut(DecNtk, cur3inpLut, 3);
    pDecLut2->nVars = 3;

    // create the right fan-in relationships
    for (int l = 0; l < 3; l++) {
        pDecLut2->vRealFis[l] = ci_sequence[l + 1];
    }

    pDecLut2->nRealFis = 3;

    return pDecLut2;  // return created and initialized If_DecLut_t.
}
int dec_l4_s3_1( If_DecSubNtk_t * DecNtk, word four_inp_lut1, word four_inp_lut2, int * vLeaves,int is_complement )
{
    int num_bits = 8;
    int size;
    int * ci_sequence;
    int ** ci_sequences = reorderArrays3( vLeaves );
    word* lut3sequence = generateBinarySequence(num_bits, &size);
    word cur3inpLut, cur3inpLutShifted;
    for (int i = 0; i < size; ++i)
    {
        int found_dec = 0;
        cur3inpLut = convertTo64Bit(lut3sequence[i], num_bits );
        for(int iter_var=0; iter_var < 4; iter_var++)
        {
            ci_sequence = ci_sequences[iter_var];
            cur3inpLutShifted = cur3inpLut;
            Abc_TtShift(&cur3inpLutShifted, 1, 4, iter_var);
            // check for correct decomposition for the first four input LUT
            for (int j = 0; j < 16; j++)
            {
                found_dec = apply_operations_3_2(four_inp_lut1, cur3inpLutShifted, (word) j, iter_var);
                if (found_dec)
                {
                    // check for correct decomposition for the second four input LUT
                    for (int k = 0; k < 16; k++)
                    {
                        found_dec = apply_operations_3_2(four_inp_lut2, cur3inpLutShifted, (word) k, iter_var);
                        if (found_dec)
                        {
                            // handle the shared three input LUT
                            If_DecLut_t *pDecLut2 = create_and_init_dec_lut2( DecNtk, cur3inpLut, ci_sequence );

                            // handle the L two input LUT
                            create_and_init_dec_lut02( DecNtk, j, is_complement, ci_sequence, pDecLut2 );

                            // handle the R two input LUT
                            create_and_init_dec_lut12( DecNtk, k, is_complement, ci_sequence, pDecLut2 );

                            // DEBUG
                           /* printf("Bit1 shift 0: \n");
                            word bit1 = s_Truths6[3];
                            word bit2 = s_Truths6[2] ;
                            word bit3 = s_Truths6[1] ;
                            word bit4 = s_Truths6[0] ;
                            printBits(bit1);
                            Abc_TtShift(&bit1, 1, 4, 2);
                            printBits(bit1);*/
                            // printf("Legal Decomposition\n");
                            /*printf("iter_var: %i\n", iter_var);
                            printf("Four input LUTs: \n");
                            printBits(four_inp_lut1);
                            printBits(four_inp_lut2);
                            printf("Three input LUT: \n");
                            printBits(cur3inpLut);
                            printf("Three input LUT shifted: \n");
                            printBits(cur3inpLutShifted);
                            printf("Two input LUTs: \n");
                            printBits(convertTo64Bit((word) j, 4));
                            printBits(convertTo64Bit((word) k, 4));
                            printf("Ci_sequence: ");
                            for (int l = 0; l < 4; l++) {
                                printf("%i ", ci_sequence[l]); // Free memory for each reordered array
                            }
                            printf("\n");*/

                            // free dynamically allocated memory of sequences
                            freeSequences1( lut3sequence, ci_sequences );
                            return 1;
                        }
                    }
                }
            }
        }
    }
    // free dynamically allocated memory of sequences
    freeSequences1( lut3sequence, ci_sequences );
    return 0;
}
/**Function*************************************************************

  Synopsis    []

  Description [Traversal functions to find shared nodes in AIG.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
int dec_l4_s3_2( If_DecSubNtk_t * DecNtk, word four_inp_lut1, word four_inp_lut2, int * vLeaves, int is_complement )
{
    int num_bits = 8;
    int size;
    int * ci_sequence;
    int ** ci_sequences = reorderArrays3( vLeaves );
    word* lut3sequence = generateBinarySequence(num_bits, &size);
    word cur3inpLut0, cur3inpLutShifted0, cur3inpLut1;
    for ( int i = 0; i < size; ++i )
    {
        int found_dec = 0;
        cur3inpLut0 = convertTo64Bit(lut3sequence[i], num_bits );
        for(int iter_var=0; iter_var < 4; iter_var++)
        {
            ci_sequence = ci_sequences[iter_var];
            cur3inpLutShifted0 = cur3inpLut0;
            Abc_TtShift(&cur3inpLutShifted0, 1, 4, iter_var);
            for ( int j = 0; j < 16; j++ )
            {
                found_dec = apply_operations_3_2(four_inp_lut1, cur3inpLutShifted0, (word) j, iter_var);
                if (found_dec)
                {
                    for ( int k = 0; k < size; ++k )
                    {
                        cur3inpLut1 = convertTo64Bit(lut3sequence[k], num_bits );
                        found_dec = apply_operations_3_3(four_inp_lut2, cur3inpLutShifted0, cur3inpLut1, iter_var);
                        if (found_dec)
                        {
                            // handle the shared three input LUT
                            If_DecLut_t *pDecLut2 = create_and_init_dec_lut2( DecNtk, cur3inpLut0, ci_sequence );
                            // handle the L two input LUT
                            create_and_init_dec_lut02( DecNtk, j, is_complement, ci_sequence, pDecLut2 );
                            // handle the R three input LUT
                            create_and_init_dec_lut13( DecNtk, cur3inpLut1, is_complement, ci_sequence, pDecLut2, vLeaves, found_dec );
                            // free dynamically allocated memory of sequences
                            freeSequences1( lut3sequence, ci_sequences );
                            return 1;
                        }
                    }
                }
            }
        }
    }
    // free memory of sequences
    freeSequences1( lut3sequence, ci_sequences );
    return 0;
}
int dec_l4_s3_3( If_DecSubNtk_t * DecNtk, word four_inp_lut1, word four_inp_lut2, int * vLeaves, int is_complement )
{
    int num_bits = 8;
    int num_bits4 = 16;
    int size, size4;
    int * ci_sequence;
    int ** ci_sequences = reorderArrays3( vLeaves );
    word* lut3sequence = generateBinarySequence(num_bits, &size);
    word* lut4sequence = generateBinarySequence(num_bits4, &size4);
    word cur3inpLut, cur3inpLutShifted, cur4inpLut;
    for ( int i = 0; i < size; ++i )
    {
        int found_dec = 0;
        cur3inpLut = convertTo64Bit(lut3sequence[i], num_bits );
        for(int iter_var=0; iter_var < 4; iter_var++)
        {
            ci_sequence = ci_sequences[iter_var];
            cur3inpLutShifted = cur3inpLut;
            Abc_TtShift(&cur3inpLutShifted, 1, 4, iter_var);
            for ( int j = 0; j < 16; j++ )
            {
                found_dec = apply_operations_3_2(four_inp_lut1, cur3inpLutShifted, (word) j, iter_var);
                if (found_dec)
                {
                    for ( int k = 0; k < size4; ++k )
                    {
                        cur4inpLut = lut4sequence[k];
                        found_dec = apply_operations_3_4(four_inp_lut2, cur3inpLutShifted, cur4inpLut, iter_var);
                        if (found_dec)
                        {
                            // handle the shared three input LUT
                            If_DecLut_t *pDecLut2 = create_and_init_dec_lut2( DecNtk, cur3inpLut, ci_sequence );
                            // handle the L two input LUT
                            create_and_init_dec_lut02( DecNtk, j, is_complement, ci_sequence, pDecLut2 );
                            // handle the R four input LUT
                            create_and_init_dec_lut14(DecNtk, cur4inpLut, is_complement, ci_sequence, pDecLut2, vLeaves, found_dec);
                            // free dynamically allocated memory of sequences
                            freeSequences2( lut3sequence, lut4sequence, ci_sequences );
                            return 1;
                        }
                    }
                }
            }
        }
    }
    // free dynamically allocated memory of sequences
    freeSequences2( lut3sequence, lut4sequence, ci_sequences );
    return 0;
}
int dec_l4_s3_4( If_DecSubNtk_t * DecNtk, word four_inp_lut1, word four_inp_lut2, int * vLeaves, int is_complement )
{
    int num_bits = 8;
    int size;
    int * ci_sequence;
    int ** ci_sequences = reorderArrays3( vLeaves );
    word* lut3sequence = generateBinarySequence(num_bits, &size);
    word cur3inpLut0, cur3inpLutShifted0, cur3inpLut1, cur3inpLut2;
    for ( int i = 0; i < size; ++i )
    {
        int found_dec0, found_dec1 = 0;
        cur3inpLut0 = convertTo64Bit(lut3sequence[i], num_bits );
        for(int iter_var=0; iter_var < 4; iter_var++)
        {
            ci_sequence = ci_sequences[iter_var];
            cur3inpLutShifted0 = cur3inpLut0;
            Abc_TtShift(&cur3inpLutShifted0, 1, 4, iter_var);
            for ( int j = 0; j < size; j++ )
            {
                cur3inpLut1 = convertTo64Bit(lut3sequence[j], num_bits );
                found_dec0 = apply_operations_3_3(four_inp_lut1, cur3inpLutShifted0, cur3inpLut1, iter_var);
                if (found_dec0)
                {
                    for ( int k = 0; k < size; ++k )
                    {
                        cur3inpLut2 = convertTo64Bit(lut3sequence[k], num_bits );
                        found_dec1 = apply_operations_3_3(four_inp_lut2, cur3inpLutShifted0, cur3inpLut2, iter_var);
                        if (found_dec1)
                        {
                            // handle the shared three input LUT
                            If_DecLut_t *pDecLut2 = create_and_init_dec_lut2( DecNtk, cur3inpLut0, ci_sequence );
                            // handle the L three input LUT
                            create_and_init_dec_lut03( DecNtk, cur3inpLut1, is_complement, ci_sequence, pDecLut2, vLeaves, found_dec0 );
                            // handle the R three input LUT
                            create_and_init_dec_lut13( DecNtk, cur3inpLut2, is_complement, ci_sequence, pDecLut2, vLeaves, found_dec1 );
                            // free dynamically allocated memory of sequences
                            freeSequences1( lut3sequence, ci_sequences );
                            return 1;
                        }
                    }
                }
            }
        }
    }
    // free memory of sequences
    freeSequences1( lut3sequence, ci_sequences );
    return 0;
}
// these solutions are already found with a shared 3 input LUT table. This is not optimal, but for the solutions it is not relevant atm.
int dec_l4_s2_1( If_DecSubNtk_t * DecNtk, word four_inp_lut1, word four_inp_lut2, int * vLeaves, int is_complement )
{
    int num_bits = 8;
    int size;
    int * ci_sequence;
    int ** ci_sequences = reorderArrays3( vLeaves );
    word* lut3sequence = generateBinarySequence(num_bits, &size);
    word cur2inpLut, cur2inpLutShifted, cur3inpLut0, cur3inpLut1;
    // 2 LUT
    for ( int i = 0; i < 16; ++i )
    {
        int found_dec = 0;
        // CHANGE for lut2sequence
        cur2inpLut = convertTo64Bit(lut3sequence[i], num_bits );
        // Iterations of 2LUT: 6 possibilities 0011 0101 0110 1001 1010 1100
        for(int iter_var=0; iter_var < 6; iter_var++)
        {
            // CHANGE sequence
            ci_sequence = ci_sequences[iter_var];
            cur2inpLutShifted = cur2inpLut;
            Abc_TtShift(&cur2inpLutShifted, 1, 4, iter_var);
            for ( int j = 0; j < size; j++ )
            {
                cur3inpLut0 = convertTo64Bit(lut3sequence[j], num_bits );
                found_dec = apply_operations_2_3(four_inp_lut1, cur2inpLutShifted, cur3inpLut0, iter_var);
                if (found_dec)
                {
                    for ( int k = 0; k < size; ++k )
                    {
                        cur3inpLut1 = convertTo64Bit(lut3sequence[k], num_bits );
                        found_dec = apply_operations_2_3(four_inp_lut2, cur2inpLutShifted, cur3inpLut1, iter_var);
                        if (found_dec)
                        {
                            // handle the shared three input LUT
                            If_DecLut_t *pDecLut2 = create_and_init_dec_lut2( DecNtk, cur3inpLut0, ci_sequence );
                            // handle the L two input LUT
                            create_and_init_dec_lut02( DecNtk, j, is_complement, ci_sequence, pDecLut2 );
                            // handle the R three input LUT
                            create_and_init_dec_lut13( DecNtk, cur3inpLut1, is_complement, ci_sequence, pDecLut2, vLeaves, found_dec );
                            // free dynamically allocated memory of sequences
                            freeSequences1( lut3sequence, ci_sequences );
                            return 1;
                        }
                    }
                }
            }
        }
    }
    // free memory of sequences
    freeSequences1( lut3sequence, ci_sequences );
    return 0;
}
/// Three shared leaves
int dec_l3_s3_1( If_DecSubNtk_t * DecNtk, word four_inp_lut1, word four_inp_lut2, int * vLeaves,int is_complement, int ind1, int ind2 )
{
    int num_bits = 8;
    int size;
    int * ci_sequence;
    int ** ci_sequences = reorderArrays3( vLeaves );
    word* lut3sequence = generateBinarySequence(num_bits, &size);
    word cur3inpLut, cur3inpLutShifted;
    for (int i = 0; i < size; ++i)
    {
        int found_dec = 0;
        cur3inpLut = convertTo64Bit(lut3sequence[i], num_bits );
        for(int iter_var=0; iter_var < 4; iter_var++)
        {
            ci_sequence = ci_sequences[iter_var];
            cur3inpLutShifted = cur3inpLut;
            Abc_TtShift(&cur3inpLutShifted, 1, 4, iter_var);
            // check for correct decomposition for the first four input LUT
            for (int j = 0; j < 16; j++)
            {
                found_dec = apply_operations_3_2(four_inp_lut1, cur3inpLutShifted, (word) j, iter_var);
                if (found_dec)
                {
                    // check for correct decomposition for the second four input LUT
                    for (int k = 0; k < 16; k++)
                    {
                        found_dec = apply_operations_3_2(four_inp_lut2, cur3inpLutShifted, (word) k, iter_var);
                        if (found_dec)
                        {
                            // handle the shared three input LUT
                            If_DecLut_t *pDecLut2 = create_and_init_dec_lut2( DecNtk, cur3inpLut, ci_sequence );

                            // handle the L two input LUT
                            create_and_init_dec_lut02( DecNtk, j, is_complement, ci_sequence, pDecLut2 );

                            // handle the R two input LUT
                            create_and_init_dec_lut12( DecNtk, k, is_complement, ci_sequence, pDecLut2 );

                            // DEBUG
                            /* printf("Bit1 shift 0: \n");
                             word bit1 = s_Truths6[3];
                             word bit2 = s_Truths6[2] ;
                             word bit3 = s_Truths6[1] ;
                             word bit4 = s_Truths6[0] ;
                             printBits(bit1);
                             Abc_TtShift(&bit1, 1, 4, 2);
                             printBits(bit1);*/
                            // printf("Legal Decomposition\n");
                            /*printf("iter_var: %i\n", iter_var);
                            printf("Four input LUTs: \n");
                            printBits(four_inp_lut1);
                            printBits(four_inp_lut2);
                            printf("Three input LUT: \n");
                            printBits(cur3inpLut);
                            printf("Three input LUT shifted: \n");
                            printBits(cur3inpLutShifted);
                            printf("Two input LUTs: \n");
                            printBits(convertTo64Bit((word) j, 4));
                            printBits(convertTo64Bit((word) k, 4));
                            printf("Ci_sequence: ");
                            for (int l = 0; l < 4; l++) {
                                printf("%i ", ci_sequence[l]); // Free memory for each reordered array
                            }
                            printf("\n");*/

                            // free dynamically allocated memory of sequences
                            freeSequences1( lut3sequence, ci_sequences );
                            return 1;
                        }
                    }
                }
            }
        }
    }
    // free dynamically allocated memory of sequences
    freeSequences1( lut3sequence, ci_sequences );
    return 0;
}

If_DecSubNtk_t* createAndInitializeSubNtk( int id, int found_id ) {
    If_DecSubNtk_t* DecNtk = (If_DecSubNtk_t*) malloc(sizeof(If_DecSubNtk_t));

    if (DecNtk == NULL) {
        printf("Memory allocation failed\n");
        return NULL;
    }

    memset(DecNtk, 0, sizeof(If_DecSubNtk_t));
    DecNtk->vLut = Vec_PtrAlloc(10);
    DecNtk->Id = id;
    DecNtk->SharedId = found_id;

    return DecNtk;
}
/**Function*************************************************************

  Synopsis    []

  Description [Traversal functions to find shared nodes in AIG.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
If_DecSubNtk_t * Abc_DecRecordToHopCros4( Abc_Ntk_t * pNtk, Hop_Man_t * pMan, If_Man_t * pIfMan, If_Cut_t * pCutBest, If_Obj_t * pIfObj, If_Obj_t * foundNode, Vec_Int_t * vCover, Vec_Int_t * vUsedNodes )
{
    Abc_Obj_t * pNode;
    If_Cut_t * pCutBest1;
    If_DecSubNtk_t * DecNtk;
    Abc_Obj_t * pNodeNew;
    int Id, SharedId;
    // Allocate memory for the Decomposed network
    DecNtk = createAndInitializeSubNtk( pIfObj->Id, foundNode->Id );

    int i;
    // if one of the nodes is already part of a subnetwork return NULL
    // the first condition ( entry == pIfObj->Id ) should always be false; if a decomposition for a network exists it is an RLut of a decomposition
    int entry;
    Vec_IntForEachEntry(vUsedNodes, entry, i) {
        if( entry == pIfObj->Id || entry == foundNode->Id )
        {
            freeDecNtk(DecNtk);
            DecNtk = NULL;
            return DecNtk;
        }
    }
    // this check should be moved somewhere inside the traverse_rec function; this should also be tried to avoid since combinations are checked multiple times
    // when commented out check if issues arise; if actually solutions are found with foundNode->Id < pIfObj->Id: something is wrong
    /*if ( foundNode->Id < pIfObj->Id   )
    {
        // printf("Id to small\n");
        freeDecNtk(DecNtk);
        DecNtk = NULL;
        return DecNtk;
    }*/
    pCutBest1 = If_ObjCutBest( foundNode );
    // Rotate pins of the other LUT as well, so they are aligned ( modify for the case of <4 shared leaves )
    /*If_Obj_t * pLeafD;
    float PinDelays[IF_MAX_LUTSIZE];
    int truthId;
    If_CutForEachLeaf( pIfMan, pCutBest1, pLeafD, i )
    {
        PinDelays[i] = If_ObjCutBest(pLeafD)->Delay;
        printf("Delay: %f\n", If_ObjCutBest(pLeafD)->Delay);
    }*/
    If_CutRotatePins( pIfMan, pCutBest1 );

    // check the leaf order
    int k;
    int LeafId;
    If_Obj_t * pLeaf;
    Abc_Obj_t * NameNode;
    int vLeaves0[4], vLeaves1[4];
    If_CutForEachLeaf( pIfMan, pCutBest, pLeaf, k )
    {
        pNodeNew = (Abc_Obj_t *)If_ObjCopy( pLeaf );
        if( pNodeNew )
        {
            LeafId = pNodeNew->Id;
            // printf( "Leaf NodeId: %i\n", pNodeNew->Id );
        }
        else
            LeafId = pLeaf->Id;
        // printf( "Leaf IfId: %i\n", pLeaf->Id );
        vLeaves0[3-k] = LeafId;
        NameNode = Abc_NtkObj( pNtk, pLeaf->Id );
        // printf("Node %i has name %s\n", NameNode->Id, Abc_ObjName(NameNode));
    }
    If_CutForEachLeaf( pIfMan, pCutBest1, pLeaf, k )
    {
        pNodeNew = (Abc_Obj_t *)If_ObjCopy( pLeaf );
        if( pNodeNew )
        {
            LeafId = pNodeNew->Id;
        }
        else
            LeafId = pLeaf->Id;
        vLeaves1[3-k] = LeafId;
        assert( vLeaves1[3-k] == vLeaves0[3-k] );
    }
    /*Abc_Obj_t * pCi;
    int x;
    Abc_NtkForEachCi(pNtk, pCi, x )
    {
        printf("Ci Ids: %i\n", pCi->Id);
    }*/

    // FIX the order of leaves
    /*for(int o = 0; o < 3; o++) {
        if(vLeaves0[i] < vLeaves0[o+1]) {  // If current element is greater than next
            // printf("No ascending order\n");
            freeDecNtk(DecNtk);
            DecNtk = NULL;
            return DecNtk;
        }
    }*/

    // maybe issue when Rotate pins really rotates something and the order gets mixed up
    /*for( int r = 0; r < 4; ++r )
    {
        printf("Leaf Ids: %i\n", vLeaves0[r]);
    }*/

    // printf( "The nodes %i and %i have shared Leaves\n", pIfObj->Id, foundNode->Id );

    // get the truth tables
    word * pTruth1 = If_CutTruthW(pIfMan, pCutBest1);
    word pTruth1_copy = *pTruth1;

    word * pTruth0 = If_CutTruthW(pIfMan, pCutBest);
    word pTruth0_copy = *pTruth0;
    // both truth tables need to be in their non-complemented form
    int is_complement = 0;
    if( If_CutTruthIsCompl( pCutBest ) )
    {
        //pCutBest->iCutFunc = 0;
        pTruth0_copy = ~pTruth0_copy;
        is_complement = 1;
    }
    if(  If_CutTruthIsCompl( pCutBest1 ) )
    {
        // pCutBest1->iCutFunc = 0;
        pTruth1_copy = ~ pTruth1_copy;
        if ( is_complement == 1 )
        {
            is_complement = 3;
        }
        else
        {
            is_complement = 2;
        }
    }
    // decomposition: l = shared leaves, s = number of inputs for the shared LUT
    // Iterate over decompositions; Start with the best and return if found.
    int (*decompose[4])(If_DecSubNtk_t *, word, word, int *, int) = {&dec_l4_s3_1, &dec_l4_s3_2, &dec_l4_s3_3, &dec_l4_s3_4 };
    int found = 0;  // Decomposition found indicator
    for( int dec_func = 0; dec_func < 4 && !found; dec_func++ )  // Loop through each decomposition function
    {
        // perform LUT-decomposition and return the LUT-structure
        if( (*decompose[dec_func])(DecNtk, pTruth0_copy, pTruth1_copy, vLeaves0, is_complement) )
        {
            // If decomposition is found, set the indicator to true and break the inner loop
            // printf("Decomposition with dec%i on root nodes: %i and %i\n", dec_func + 1, pIfObj->Id, foundNode->Id );
            found = 1;
            break;
        }
    }
    if ( found == 0 )
    {
        freeSubnetwork(DecNtk);
        // printf( "No decomposition found\n" );
        DecNtk = NULL;
    }
    assert( !pIfMan->pPars->fUseTtPerm );
    return DecNtk;
}
// Function to check if a number exists in the array
int is_exist(int num, int arr[], int size) {
    for(int i = 0; i < size; ++i) {
        if(arr[i] == num)
            return 1;
    }
    return 0;
}
// Function to find the position of different element
int find_diff_position(int arr1[], int arr2[], int size) {
    for(int i = 0; i < size; ++i) {
        if(!is_exist(arr1[i], arr2, size))
            return i;
    }
    return -1; // return -1 if no difference found
}
If_DecSubNtk_t * Abc_DecRecordToHopCros3( Abc_Ntk_t * pNtk, Hop_Man_t * pMan, If_Man_t * pIfMan, If_Cut_t * pCutBest, If_Obj_t * pIfObj, If_Obj_t * foundNode, Vec_Int_t * vCover, Vec_Int_t * vUsedNodes )
{
    Abc_Obj_t * pNode;
    If_Cut_t * pCutBest1;
    If_DecSubNtk_t * DecNtk;
    Abc_Obj_t * pNodeNew;
    int Id, SharedId;
    // Allocate memory for the Decomposed network
    DecNtk = createAndInitializeSubNtk( pIfObj->Id, foundNode->Id );

    printf("Three shared Leaves\n");

    int i;
    // if one of the nodes is already part of a subnetwork return NULL
    // the first condition ( entry == pIfObj->Id ) should always be false; if a decomposition for a network exists it is an RLut of a decomposition
    int entry;
    Vec_IntForEachEntry(vUsedNodes, entry, i) {
        if( entry == pIfObj->Id || entry == foundNode->Id )
        {
            freeDecNtk(DecNtk);
            DecNtk = NULL;
            return DecNtk;
        }
    }
    pCutBest1 = If_ObjCutBest( foundNode );
    // Rotate pins of the other LUT as well, so they are aligned ( modify for the case of <4 shared leaves )
    /*If_Obj_t * pLeafD;
    float PinDelays[IF_MAX_LUTSIZE];
    int truthId;
    If_CutForEachLeaf( pIfMan, pCutBest1, pLeafD, i )
    {
        PinDelays[i] = If_ObjCutBest(pLeafD)->Delay;
        printf("Delay: %f\n", If_ObjCutBest(pLeafD)->Delay);
    }*/
    If_CutRotatePins( pIfMan, pCutBest1 );

    // check the leaf order
    int k;
    int LeafId;
    If_Obj_t * pLeaf;
    Abc_Obj_t * NameNode;
    int vLeaves0[4], vLeaves1[4];
    printf("Leaves of first cut\n");
    If_CutForEachLeaf( pIfMan, pCutBest, pLeaf, k )
    {
        pNodeNew = (Abc_Obj_t *)If_ObjCopy( pLeaf );
        if( pNodeNew )
        {
            LeafId = pNodeNew->Id;
            printf( "Leaf NodeId: %i\n", pNodeNew->Id );
        }
        else
            LeafId = pLeaf->Id;
        printf( "Leaf IfId: %i\n", pLeaf->Id );
        vLeaves0[3-k] = LeafId;
        NameNode = Abc_NtkObj( pNtk, pLeaf->Id );
        // printf("Node %i has name %s\n", NameNode->Id, Abc_ObjName(NameNode));
    }
    printf("Leaves of second cut\n");
    If_CutForEachLeaf( pIfMan, pCutBest1, pLeaf, k )
    {
        pNodeNew = (Abc_Obj_t *)If_ObjCopy( pLeaf );
        if( pNodeNew )
        {
            LeafId = pNodeNew->Id;
            printf( "Leaf NodeId: %i\n", pNodeNew->Id );
        }
        else
            LeafId = pLeaf->Id;
        printf( "Leaf IfId: %i\n", pLeaf->Id );
        vLeaves1[3-k] = LeafId;
        // assert( vLeaves1[3-k] == vLeaves0[3-k] );
    }
    int index1 = find_diff_position(vLeaves0, vLeaves1, 4);
    int index2 = find_diff_position(vLeaves1, vLeaves0, 4);
    if(index1 != -1) {
        printf("The different element in arr1 is at position %d\n", index1);
    } else {
        printf("No different elements found in arr1\n");
    }

    if(index2 != -1) {
        printf("The different element in arr2 is at position %d\n", index2);
    } else {
        printf("No different elements found in arr2\n");
    }

    // get the truth tables
    word * pTruth1 = If_CutTruthW(pIfMan, pCutBest1);
    word pTruth1_copy = *pTruth1;

    word * pTruth0 = If_CutTruthW(pIfMan, pCutBest);
    word pTruth0_copy = *pTruth0;
    // both truth tables need to be in their non-complemented form
    int is_complement = 0;
    if( If_CutTruthIsCompl( pCutBest ) )
    {
        //pCutBest->iCutFunc = 0;
        pTruth0_copy = ~pTruth0_copy;
        is_complement = 1;
    }
    if(  If_CutTruthIsCompl( pCutBest1 ) )
    {
        // pCutBest1->iCutFunc = 0;
        pTruth1_copy = ~ pTruth1_copy;
        if ( is_complement == 1 )
        {
            is_complement = 3;
        }
        else
        {
            is_complement = 2;
        }
    }
    // decomposition: l = shared leaves, s = number of inputs for the shared LUT
    // Iterate over decompositions; Start with the best and return if found.
    int (*decompose[1])(If_DecSubNtk_t *, word, word, int *, int, int, int) = {&dec_l3_s3_1 };
    int found = 0;  // Decomposition found indicator
    for( int dec_func = 0; dec_func < 1 && !found; dec_func++ )  // Loop through each decomposition function
    {
        // perform LUT-decomposition and return the LUT-structure
        if( (*decompose[dec_func])(DecNtk, pTruth0_copy, pTruth1_copy, vLeaves0, is_complement, index1, index2) )
        {
            // If decomposition is found, set the indicator to true and break the inner loop
            // printf("Decomposition with dec%i on root nodes: %i and %i\n", dec_func + 1, pIfObj->Id, foundNode->Id );
            found = 1;
            break;
        }
    }
    if ( found == 0 )
    {
        freeSubnetwork(DecNtk);
        // printf( "No decomposition found\n" );
        DecNtk = NULL;
    }
    assert( !pIfMan->pPars->fUseTtPerm );
    return DecNtk;
}

/**Function*************************************************************

  Synopsis    []

  Description [Unused function: Should restore the right Obj order in the Vec_Prt holding vObjs in the Ntk]

  SideEffects [Also the fan-in relationships have to be modified for this to work.]

  SeeAlso     []

***********************************************************************/
void restore_order( Abc_Ntk_t * pNtkNew, If_Obj_t * pIfObj )
{
    void* lastElement = Vec_PtrEntry(pNtkNew->vObjs, Vec_PtrSize(pNtkNew->vObjs) - 1);

    // Shift elements in vector from desired index to end
    for(int idx = Vec_PtrSize(pNtkNew->vObjs) - 1; idx > pIfObj->Id; --idx) {
        Vec_PtrWriteEntry(pNtkNew->vObjs, idx, Vec_PtrEntry(pNtkNew->vObjs, idx - 1));
    }

    // Write the temporary element at the desired position
    Vec_PtrWriteEntry(pNtkNew->vObjs, pIfObj->Id, lastElement);
}

/**********************************************************************************************************************/
/**********************************************************************************************************************/
/**********************************************************************************************************************/


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
/**********************************************************************************************************************/
/**********************************************************************************************************************/
/**********************************************************************************************************************/


#endif //ABC_IFCROSDEC_H
