moris::Mat<moris::uint> tAnc1DInds(1,50,UINT_MAX);
moris::Mat<moris::uint> tAnc1DRank(1,50,UINT_MAX);


// Edges created on stk edge
tAnc1DInds(0,0)  = aParentEntities(1)(0,0);  // XTK edges created on parent element local edge 1
tAnc1DInds(0,1)  = aParentEntities(1)(0,1);  // XTK edges created on parent element local edge 2
tAnc1DInds(0,2)  = aParentEntities(1)(0,2);  // XTK edges created on parent element local edge 3
tAnc1DInds(0,3)  = aParentEntities(1)(0,3);  // XTK edges created on parent element local edge 4
tAnc1DInds(0,4)  = aParentEntities(1)(0,4);  // XTK edges created on parent element local edge 5
tAnc1DInds(0,5)  = aParentEntities(1)(0,5);  // XTK edges created on parent element local edge 6
tAnc1DInds(0,6)  = aParentEntities(1)(0,6);  // XTK edges created on parent element local edge 7
tAnc1DInds(0,7)  = aParentEntities(1)(0,7);  // XTK edges created on parent element local edge 8
tAnc1DInds(0,8)  = aParentEntities(1)(0,8);  // XTK edges created on parent element local edge 9
tAnc1DInds(0,9)  = aParentEntities(1)(0,9);  // XTK edges created on parent element local edge 10
tAnc1DInds(0,10) = aParentEntities(1)(0,10); // XTK edges created on parent element local edge 11
tAnc1DInds(0,11) = aParentEntities(1)(0,11); // XTK edges created on parent element local edge 12

// Edges created within stk element
tAnc1DInds(0,12)  = mParentEID;  // XTK edges created on parent element local edge 13
tAnc1DInds(0,13)  = mParentEID;  // XTK edges created on parent element local edge 14
tAnc1DInds(0,14)  = mParentEID;  // XTK edges created on parent element local edge 15
tAnc1DInds(0,15)  = mParentEID;  // XTK edges created on parent element local edge 16
tAnc1DInds(0,16)  = mParentEID;  // XTK edges created on parent element local edge 17
tAnc1DInds(0,17)  = mParentEID;  // XTK edges created on parent element local edge 18
tAnc1DInds(0,18)  = mParentEID;  // XTK edges created on parent element local edge 19
tAnc1DInds(0,19)  = mParentEID;  // XTK edges created on parent element local edge 20
tAnc1DInds(0,44)  = mParentEID;  // XTK edges created on parent element local edge 45
tAnc1DInds(0,45)  = mParentEID;  // XTK edges created on parent element local edge 46
tAnc1DInds(0,46)  = mParentEID;  // XTK edges created on parent element local edge 47
tAnc1DInds(0,47)  = mParentEID;  // XTK edges created on parent element local edge 48
tAnc1DInds(0,48)  = mParentEID;  // XTK edges created on parent element local edge 49
tAnc1DInds(0,49)  = mParentEID;  // XTK edges created on parent element local edge 50

//Edges created on stk face
tAnc1DInds(0,20) = aParentEntities(2)(0,2); // XTK edges created on parent element local edge 21
tAnc1DInds(0,21) = aParentEntities(2)(0,2); // XTK edges created on parent element local edge 22
tAnc1DInds(0,22) = aParentEntities(2)(0,2); // XTK edges created on parent element local edge 23
tAnc1DInds(0,23) = aParentEntities(2)(0,2); // XTK edges created on parent element local edge 24

tAnc1DInds(0,24) = aParentEntities(2)(0,1); // XTK edges created on parent element local edge 25
tAnc1DInds(0,25) = aParentEntities(2)(0,1); // XTK edges created on parent element local edge 26
tAnc1DInds(0,26) = aParentEntities(2)(0,1); // XTK edges created on parent element local edge 27
tAnc1DInds(0,27) = aParentEntities(2)(0,1); // XTK edges created on parent element local edge 28

tAnc1DInds(0,28) = aParentEntities(2)(0,0); // XTK edges created on parent element local edge 29
tAnc1DInds(0,29) = aParentEntities(2)(0,0); // XTK edges created on parent element local edge 30
tAnc1DInds(0,30) = aParentEntities(2)(0,0); // XTK edges created on parent element local edge 31
tAnc1DInds(0,31) = aParentEntities(2)(0,0); // XTK edges created on parent element local edge 32

tAnc1DInds(0,32) = aParentEntities(2)(0,3); // XTK edges created on parent element local edge 33
tAnc1DInds(0,33) = aParentEntities(2)(0,3); // XTK edges created on parent element local edge 34
tAnc1DInds(0,34) = aParentEntities(2)(0,3); // XTK edges created on parent element local edge 35
tAnc1DInds(0,35) = aParentEntities(2)(0,3); // XTK edges created on parent element local edge 36

tAnc1DInds(0,36) = aParentEntities(2)(0,5);  // XTK edges created on parent element local edge 37
tAnc1DInds(0,37) = aParentEntities(2)(0,5);  // XTK edges created on parent element local edge 38
tAnc1DInds(0,38) = aParentEntities(2)(0,5);  // XTK edges created on parent element local edge 39
tAnc1DInds(0,39) = aParentEntities(2)(0,5); // XTK edges created on parent element local edge 40

tAnc1DInds(0,40) = aParentEntities(2)(0,4); // XTK edges created on parent element local edge 41
tAnc1DInds(0,41) = aParentEntities(2)(0,4); // XTK edges created on parent element local edge 42
tAnc1DInds(0,42) = aParentEntities(2)(0,4); // XTK edges created on parent element local edge 43
tAnc1DInds(0,43) = aParentEntities(2)(0,4); // XTK edges created on parent element local edge 44




// Edges created on stk edge
tAnc1DRank(0,0)  = 1; // XTK edges created on parent element local edge 1
tAnc1DRank(0,1)  = 1; // XTK edges created on parent element local edge 2
tAnc1DRank(0,2)  = 1; // XTK edges created on parent element local edge 3
tAnc1DRank(0,3)  = 1; // XTK edges created on parent element local edge 4
tAnc1DRank(0,4)  = 1; // XTK edges created on parent element local edge 5
tAnc1DRank(0,5)  = 1; // XTK edges created on parent element local edge 6
tAnc1DRank(0,6)  = 1; // XTK edges created on parent element local edge 7
tAnc1DRank(0,7)  = 1; // XTK edges created on parent element local edge 8
tAnc1DRank(0,8)  = 1; // XTK edges created on parent element local edge 9
tAnc1DRank(0,9)  = 1; // XTK edges created on parent element local edge 10
tAnc1DRank(0,10) = 1; // XTK edges created on parent element local edge 11
tAnc1DRank(0,11) = 1; // XTK edges created on parent element local edge 12

// Edges created within stk element
tAnc1DRank(0,12)  = 3; // XTK edges created on parent element local edge 13
tAnc1DRank(0,13)  = 3; // XTK edges created on parent element local edge 14
tAnc1DRank(0,14)  = 3; // XTK edges created on parent element local edge 15
tAnc1DRank(0,15)  = 3; // XTK edges created on parent element local edge 16
tAnc1DRank(0,16)  = 3; // XTK edges created on parent element local edge 17
tAnc1DRank(0,17)  = 3; // XTK edges created on parent element local edge 18
tAnc1DRank(0,18)  = 3; // XTK edges created on parent element local edge 19
tAnc1DRank(0,19)  = 3; // XTK edges created on parent element local edge 20
tAnc1DRank(0,44)  = 3; // XTK edges created on parent element local edge 45
tAnc1DRank(0,45)  = 3; // XTK edges created on parent element local edge 46
tAnc1DRank(0,46)  = 3; // XTK edges created on parent element local edge 47
tAnc1DRank(0,47)  = 3; // XTK edges created on parent element local edge 48
tAnc1DRank(0,48)  = 3; // XTK edges created on parent element local edge 49
tAnc1DRank(0,49)  = 3; // XTK edges created on parent element local edge 50

//Edges created on stk face
tAnc1DRank(0,20) = 2; // XTK edges created on parent element local edge 21
tAnc1DRank(0,21) = 2; // XTK edges created on parent element local edge 22
tAnc1DRank(0,22) = 2; // XTK edges created on parent element local edge 23
tAnc1DRank(0,23) = 2; // XTK edges created on parent element local edge 24

tAnc1DRank(0,24) = 2; // XTK edges created on parent element local edge 25
tAnc1DRank(0,25) = 2; // XTK edges created on parent element local edge 26
tAnc1DRank(0,26) = 2; // XTK edges created on parent element local edge 27
tAnc1DRank(0,27) = 2; // XTK edges created on parent element local edge 28

tAnc1DRank(0,28) = 2; // XTK edges created on parent element local edge 29
tAnc1DRank(0,29) = 2; // XTK edges created on parent element local edge 30
tAnc1DRank(0,30) = 2; // XTK edges created on parent element local edge 31
tAnc1DRank(0,31) = 2; // XTK edges created on parent element local edge 32

tAnc1DRank(0,32) = 2; // XTK edges created on parent element local edge 33
tAnc1DRank(0,33) = 2; // XTK edges created on parent element local edge 34
tAnc1DRank(0,34) = 2; // XTK edges created on parent element local edge 35
tAnc1DRank(0,35) = 2; // XTK edges created on parent element local edge 36

tAnc1DRank(0,36) = 2;  // XTK edges created on parent element local edge 37
tAnc1DRank(0,37) = 2;  // XTK edges created on parent element local edge 38
tAnc1DRank(0,38) = 2;  // XTK edges created on parent element local edge 39
tAnc1DRank(0,39) = 2; // XTK edges created on parent element local edge 40

tAnc1DRank(0,40) = 2; // XTK edges created on parent element local edge 41
tAnc1DRank(0,41) = 2; // XTK edges created on parent element local edge 42
tAnc1DRank(0,42) = 2; // XTK edges created on parent element local edge 43
tAnc1DRank(0,43) = 2; // XTK edges created on parent element local edge 44
