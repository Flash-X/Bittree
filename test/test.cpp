#include <gtest/gtest.h>
#include <iostream>
#include <algorithm>

#include "macros.h"
#include "Bittree_core.h"

namespace {
#if NDIM==1
    int mort_true[6] = {1,2,3,4,5,6};
    int bitid_true[6] = {2,3,4,5,6,7};
#elif NDIM==2
    int mort_true[46] = {1, 2, 15, 28, 45, 46, 3, 4, 9, 10, 16, 17, 29, 30, 18, 23, 35, 40, 5, 6, 7, 8, 11, 12, 13, 14, 31, 32, 33, 34, 19, 20, 24, 25, 36, 37, 41, 42, 21, 22, 26, 27, 38, 39, 43, 44};
    int bitid_true[46] = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 21, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 40, 41, 42, 43, 32, 33, 36, 37, 44, 45, 48, 49, 34, 35, 38, 39, 46, 47, 50, 51};
#else
    int mort_true[376] = {
        1, 2, 43, 84, 361, 362, 141, 182, 239, 296, 363, 364, 365, 366, 367, 368, 373, 374, 369, 370,
        371, 372, 375, 376, 3, 4, 13, 14, 44, 45, 85, 86, 46, 55, 95, 104, 23, 24, 33, 34,
        64, 65, 113, 114, 66, 75, 123, 132, 142, 143, 183, 184, 144, 145, 193, 194, 240, 241, 297, 298,
        242, 251, 307, 316, 146, 155, 203, 212, 164, 173, 221, 230, 260, 269, 325, 334, 278, 287, 343, 352,
        5, 6, 7, 8, 15, 16, 17, 18, 87, 88, 89, 90, 47, 48, 56, 57, 96, 97, 105, 106,
        49, 50, 58, 59, 98, 99, 107, 108, 9, 10, 11, 12, 19, 20, 21, 22, 91, 92, 93, 94,
        51, 52, 60, 61, 100, 101, 109, 110, 53, 54, 62, 63, 102, 103, 111, 112, 25, 26, 27, 28,
        35, 36, 37, 38, 115, 116, 117, 118, 67, 68, 76, 77, 124, 125, 133, 134, 69, 70, 78, 79,
        126, 127, 135, 136, 29, 30, 31, 32, 39, 40, 41, 42, 119, 120, 121, 122, 71, 72, 80, 81,
        128, 129, 137, 138, 73, 74, 82, 83, 130, 131, 139, 140, 185, 186, 187, 188, 195, 196, 197, 198,
        299, 300, 301, 302, 243, 244, 252, 253, 308, 309, 317, 318, 245, 246, 254, 255, 310, 311, 319, 320,
        189, 190, 191, 192, 199, 200, 201, 202, 303, 304, 305, 306, 247, 248, 256, 257, 312, 313, 321, 322,
        249, 250, 258, 259, 314, 315, 323, 324, 147, 148, 156, 157, 204, 205, 213, 214, 149, 150, 158, 159,
        206, 207, 215, 216, 165, 166, 174, 175, 222, 223, 231, 232, 167, 168, 176, 177, 224, 225, 233, 234,
        261, 262, 270, 271, 326, 327, 335, 336, 263, 264, 272, 273, 328, 329, 337, 338, 279, 280, 288, 289,
        344, 345, 353, 354, 281, 282, 290, 291, 346, 347, 355, 356, 151, 152, 160, 161, 208, 209, 217, 218,
        153, 154, 162, 163, 210, 211, 219, 220, 169, 170, 178, 179, 226, 227, 235, 236, 171, 172, 180, 181,
        228, 229, 237, 238, 265, 266, 274, 275, 330, 331, 339, 340, 267, 268, 276, 277, 332, 333, 341, 342,
        283, 284, 292, 293, 348, 349, 357, 358, 285, 286, 294, 295, 350, 351, 359, 360};
    int bitid_true[376] = {
        24, 25, 26, 27, 32, 33, 28, 29, 30, 31, 34, 35, 36, 37, 38, 39, 44, 45, 40, 41,
        42, 43, 46, 47, 48, 49, 50, 51, 56, 57, 64, 65, 58, 59, 66, 67, 52, 53, 54, 55,
        60, 61, 68, 69, 62, 63, 70, 71, 72, 73, 80, 81, 74, 75, 82, 83, 88, 89, 96, 97,
        90, 91, 98, 99, 76, 77, 84, 85, 78, 79, 86, 87, 92, 93, 100, 101, 94, 95, 102, 103,
        104, 105, 106, 107, 112, 113, 114, 115, 168, 169, 170, 171, 136, 137, 144, 145, 176, 177, 184, 185,
        138, 139, 146, 147, 178, 179, 186, 187, 108, 109, 110, 111, 116, 117, 118, 119, 172, 173, 174, 175,
        140, 141, 148, 149, 180, 181, 188, 189, 142, 143, 150, 151, 182, 183, 190, 191, 120, 121, 122, 123,
        128, 129, 130, 131, 192, 193, 194, 195, 152, 153, 160, 161, 200, 201, 208, 209, 154, 155, 162, 163,
        202, 203, 210, 211, 124, 125, 126, 127, 132, 133, 134, 135, 196, 197, 198, 199, 156, 157, 164, 165,
        204, 205, 212, 213, 158, 159, 166, 167, 206, 207, 214, 215, 248, 249, 250, 251, 256, 257, 258, 259,
        344, 345, 346, 347, 296, 297, 304, 305, 352, 353, 360, 361, 298, 299, 306, 307, 354, 355, 362, 363,
        252, 253, 254, 255, 260, 261, 262, 263, 348, 349, 350, 351, 300, 301, 308, 309, 356, 357, 364, 365,
        302, 303, 310, 311, 358, 359, 366, 367, 216, 217, 224, 225, 264, 265, 272, 273, 218, 219, 226, 227,
        266, 267, 274, 275, 232, 233, 240, 241, 280, 281, 288, 289, 234, 235, 242, 243, 282, 283, 290, 291,
        312, 313, 320, 321, 368, 369, 376, 377, 314, 315, 322, 323, 370, 371, 378, 379, 328, 329, 336, 337,
        384, 385, 392, 393, 330, 331, 338, 339, 386, 387, 394, 395, 220, 221, 228, 229, 268, 269, 276, 277,
        222, 223, 230, 231, 270, 271, 278, 279, 236, 237, 244, 245, 284, 285, 292, 293, 238, 239, 246, 247,
        286, 287, 294, 295, 316, 317, 324, 325, 372, 373, 380, 381, 318, 319, 326, 327, 374, 375, 382, 383,
        332, 333, 340, 341, 388, 389, 396, 397, 334, 335, 342, 343, 390, 391, 398, 399};
#endif

class BittreeUnitTest : public testing::Test {
protected:
    BittreeUnitTest(void) {
      int nbase = CONCAT_NDIM(2,*3,*4);
      int top[NDIM] = {LIST_NDIM(2,3,4)};
      bool includes[nbase];
      for(int i=0; i<nbase; i++) includes[i] = true;
      bittree_init(top, includes);
    }

    ~BittreeUnitTest(void) {
        bool updated = true;
        int count;
        bittree_block_count(&updated, &count);
        std::cout << "End of test block count=" << count << std::endl;
    }
};

TEST_F(BittreeUnitTest,RefinementTest){

    // Declare some variables
    int lev,mort,bitid,count;
    int ijk[3];
    bool updated, val;
    int xlim, ylim, zlim;

    // First round of refinement
    bittree_refine_init();

    xlim = std::max(2*K1D,1);
    ylim = std::max(3*K2D,1);
    zlim = std::max(4*K3D,1);
    updated = false;
    val = true;
    bitid=0;
    mort=0;

    for(    int k=0; k<zlim; ++k) {
      for(  int j=0; j<ylim; ++j) {
        for(int i=0; i<xlim; ++i) {
          lev = 0;
          ijk[0]=i; ijk[1]=j; ijk[2]=k;
          bittree_identify(&updated,&lev,ijk,&mort,&bitid);
          if( (i==1 || j==1 || k==1) && i<=1 && j<=1 && k<=1) {
              bittree_refine_mark(&bitid,&val);
          }
    }}}

    bittree_refine_update();
    bittree_refine_apply();

    // Check block count is correct after first refinement
    updated = false;
    bittree_block_count(&updated, &count);
    ASSERT_EQ(count, SELECT_NDIM(4,18,80) );

    // Second round of refinement
    bittree_refine_init();

    xlim = std::max(4*K1D,1);
    ylim = std::max(6*K2D,1);
    zlim = std::max(8*K3D,1);
    count = 0;
    updated = false;

    for(    int k=0; k<zlim; ++k) {
      for(  int j=0; j<ylim; ++j) {
        for(int i=0; i<xlim; ++i) {
          lev = 1;
          ijk[0]=i; ijk[1]=j; ijk[2]=k;
          bittree_identify(&updated,&lev,ijk,&mort,&bitid);
          if(lev!=1) continue;

          if( (i==3 || j==3 || k==3) && i<=3 && j<=3 && k<=3) {
              bittree_refine_mark(&bitid,&val);
          }
    }}}


    bittree_refine_update();
    bittree_refine_apply();

    // Check block count is correct after second refinement
    updated = false;
    bittree_block_count(&updated, &count);
    ASSERT_EQ(count, SELECT_NDIM(6,46,376) );

    // Check tree with bittree_identify and bittree_locate
    xlim = std::max(2*K1D,1);
    ylim = std::max(3*K2D,1);
    zlim = std::max(4*K3D,1);
    updated = false;
    count = 0;

    for(int l=0; l<3; ++l) {
      if(l>0) {
          xlim = xlim*(1+K1D);
          ylim = ylim*(1+K2D);
          zlim = zlim*(1+K3D);
      }
      for(    int k=0; k<zlim; ++k) {
        for(  int j=0; j<ylim; ++j) {
          for(int i=0; i<xlim; ++i) {
            lev = l;
            ijk[0]=i; ijk[1]=j; ijk[2]=k;
            bittree_identify(&updated,&lev,ijk,&mort,&bitid);
            mort++;
            if(lev!=l) continue;

            ASSERT_EQ( mort,  mort_true[count]  );
            ASSERT_EQ( bitid, bitid_true[count] );

            count ++;
      }}}
    }


}

// Test Bittree core functions
TEST_F(BittreeUnitTest,BittreeCore){
    int lev,mort,bitid,count,id0;
    int id_lims[2];
    int ijk[3];
    bool updated, val;

    val = bittree_initialized();
    ASSERT_EQ( val, true);

    bittree_refine_init();

    // Test get_id0
    updated = false;
    bittree_get_id0(&updated, &id0);
    ASSERT_EQ( id0, SELECT_NDIM(2,6,24) );

    // Mark a single block for refinement
    bitid = id0;
    val = true;
    bittree_refine_mark(&bitid,&val);
    bittree_refine_update();

    updated = true;
    bittree_get_id0(&updated, &id0);
    ASSERT_EQ( id0, SELECT_NDIM(2,6,24) );

    // Test level count
    updated = false;
    bittree_level_count(&updated, &count);
    ASSERT_EQ( count, 1);
    updated = true;
    bittree_level_count(&updated, &count);
    ASSERT_EQ( count, 2);

    // Test block count
    updated = false;
    bittree_block_count(&updated, &count);
    ASSERT_EQ( count, SELECT_NDIM(2,6,24));
    updated = true;
    bittree_block_count(&updated, &count);
    ASSERT_EQ( count, SELECT_NDIM(4,10,32));

    // Test leaf count
    updated = false;
    bittree_leaf_count(&updated, &count);
    ASSERT_EQ( count, SELECT_NDIM(2,6,24));
    updated = true;
    bittree_leaf_count(&updated, &count);
    ASSERT_EQ( count, SELECT_NDIM(3,9,31));
    
    // Test delta count
    bittree_delta_count(&count);
    ASSERT_EQ( count, 1);

    // Test check_refine_bit
    bitid = id0;
    bittree_check_refine_bit(&bitid,&val);
    ASSERT_EQ( val, true);
    bitid = id0+1;
    bittree_check_refine_bit(&bitid,&val);
    ASSERT_EQ( val, false);

    // Test is_parent
    bitid = id0;
    updated = false;
    bittree_is_parent( &updated, &bitid, &val);
    ASSERT_EQ( val, false);
    bitid = id0;
    updated = true;
    bittree_is_parent( &updated, &bitid, &val);
    ASSERT_EQ( val, true);

    // Test bittree_identify
    //   Not updated - Wrong level
    lev = 1;
    ijk[0]=0; ijk[1]=0; ijk[2]=0;
    updated = false;
    bittree_identify(&updated,&lev,ijk,&mort,&bitid);
    ASSERT_EQ(lev,0);
    ASSERT_EQ(mort,0);
    ASSERT_EQ(bitid,SELECT_NDIM(2,6,24));

    //    Not updated - Out of bounds
    lev = 0;
    ijk[0]=234; ijk[1]=42; ijk[2]=577;
    updated = false;
    bittree_identify(&updated,&lev,ijk,&mort,&bitid);
    ASSERT_EQ(lev,-1);
    ASSERT_EQ(mort,-1);
    ASSERT_EQ(bitid,-1);

    //    Updated
    lev = 1;
    ijk[0]=0; ijk[1]=0; ijk[2]=0;
    updated = true;
    bittree_identify(&updated,&lev,ijk,&mort,&bitid);
    ASSERT_EQ(lev,1);
    ASSERT_EQ(mort,SELECT_NDIM(1,1,1) );
    ASSERT_EQ(bitid, SELECT_NDIM(4,12,48) );

    //    Updated - out of bounds
    lev = 0;
    ijk[0]=65; ijk[1]=100; ijk[2]=700;
    updated = true;
    bittree_identify(&updated,&lev,ijk,&mort,&bitid);
    ASSERT_EQ(lev,-1);
    ASSERT_EQ(mort,-1);
    ASSERT_EQ(bitid,-1);

    // Test get_level_id_limits
    updated = false;
    lev = 0;
    bittree_level_bitid_limits(&updated, &lev, id_lims);
    ASSERT_EQ( id_lims[0], SELECT_NDIM(2,6,24));
    ASSERT_EQ( id_lims[1], SELECT_NDIM(4,12,48));

    updated = true;
    lev = 1;
    bittree_level_bitid_limits(&updated, &lev, id_lims);
    ASSERT_EQ( id_lims[0], SELECT_NDIM(4,12,48));
    ASSERT_EQ( id_lims[1], SELECT_NDIM(6,16,56));

    // Test bittree_locate
    bitid = id0 + 1;
    updated = false;
    bittree_locate(&updated,&bitid,&lev,ijk,&mort);
    ASSERT_EQ(lev,0);
    ASSERT_EQ(mort,1);
    ASSERT_EQ(ijk[0], 1);
    for(int n=1; n<NDIM; ++n) ASSERT_EQ(ijk[n], 0);

    bitid = 2112;
    updated = false;
    bittree_locate(&updated,&bitid,&lev,ijk,&mort);
    ASSERT_EQ(lev,-1);
    ASSERT_EQ(mort,-1);
    for(int n=0; n<NDIM; ++n) ASSERT_EQ(ijk[n], -1);

    bitid = id_lims[0]+1;
    updated = true;
    bittree_locate(&updated,&bitid,&lev,ijk,&mort);
    ASSERT_EQ(lev,1);
    ASSERT_EQ(mort,2);
    ASSERT_EQ(ijk[0], 1);
    for(int n=1; n<NDIM; ++n) ASSERT_EQ(ijk[n], 0);

    bitid = 4532;
    updated = true;
    bittree_locate(&updated,&bitid,&lev,ijk,&mort);
    ASSERT_EQ(lev,-1);
    ASSERT_EQ(mort,-1);
    for(int n=0; n<NDIM; ++n) ASSERT_EQ(ijk[n], -1);

    // Test get_bitid_list
    updated = false;
    int mmin = 1;
    int mmax = SELECT_NDIM(2,4,19);
    int bitid_list[mmax-mmin];
#if NDIM==1
    int true_list[2] = {2,3};
#elif NDIM==2
    int true_list[6] = {6,7,8,9,10,11};
#else
    int true_list[24] = {24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47};
#endif
    bittree_get_bitid_list(&updated, &mmin, &mmax, bitid_list);
    for( int i=mmin; i<mmax; ++i) {
        ASSERT_EQ( bitid_list[i-mmin], true_list[i] ); 
    }

    updated = true;
    mmin = 1;
    mmax = SELECT_NDIM(3,7,28);
    int bitid_list_2[mmax-mmin];
#if NDIM==1
    int true_list_2[4] = {2,4,5,3};
#elif NDIM==2
    int true_list_2[10] = {6,12,13,14,15,7,8,9,10,11};
#else
    int true_list_2[32] = {24,48,49,50,51,52,53,54,55,25,26,27,28,29,30,31,32,33,34,35,36,
                           37,38,39,40,41,42,43,44,45,46,47};
#endif
    bittree_get_bitid_list(&updated, &mmin, &mmax, bitid_list_2);
    for( int i=mmin; i<mmax; ++i) {
        ASSERT_EQ( bitid_list_2[i-mmin], true_list_2[i] ); 
    }


    // Test refine_reduce
    // TODO

    // Test print_2d
    int dtype = 0;
    bittree_print_2d(&dtype);

    bittree_refine_apply();

    // Test delta count outside of refine
    bittree_delta_count(&count);
    ASSERT_EQ( count, 0);

    // Test check_refine_bit outside of refine
    bitid = id0;
    bittree_check_refine_bit(&bitid,&val);
    ASSERT_EQ( val, false);
}

}
