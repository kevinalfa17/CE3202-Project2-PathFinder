#include <gtest/gtest.h>
#include "config.h"
#include "Matrix/Matrix.hpp"
#include "PathFinder/IndexMap.h"

using namespace std;
using namespace anpi;

template class anpi::Matrix<double>;
template class anpi::Matrix<float>;

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

//Constructor
TEST(MatrixConstructors, Default) {
    
    anpi::Matrix<float> a;
    ASSERT_EQ(a.rows() , 0) << "Wrong rows";
    ASSERT_EQ(a.cols() , 0) << "Wrong cols";
    
}

TEST(MatrixConstructors, Unitialiazed) {
    
    anpi::Matrix<float> a(2,3,anpi::Matrix<float>::DoNotInitialize);
    ASSERT_EQ(a.rows() , 2) << "Wrong rows";
    ASSERT_EQ(a.cols() , 3) << "Wrong cols";
    
}

TEST(MatrixConstructors, DefaultInitialiazed1) {
    
    anpi::Matrix<float> a(3,2);
    ASSERT_EQ(a.rows() , 3) << "Wrong rows";
    ASSERT_EQ(a.cols() , 2) << "Wrong cols";
    ASSERT_EQ(a(0,0) , 0.f) << "Wrong value";
    
}

TEST(MatrixConstructors, DefaultInitialiazed2) {
    
    anpi::Matrix<double> a(3,2,4.);
    ASSERT_EQ(a.rows() , 3) << "Wrong rows";
    ASSERT_EQ(a.cols() , 2) << "Wrong cols";
    ASSERT_EQ(a(0,0) , 4.) << "Wrong value";
    
}

TEST(MatrixConstructors, UnitialiazedPadded) {
    
    anpi::Matrix<float> a(2,3,anpi::Matrix<float>::Padded);
    ASSERT_EQ(a.rows() , 2) << "Wrong rows";
    ASSERT_EQ(a.cols() , 3) << "Wrong cols";
    
}

TEST(MatrixConstructors, InitialiazedPadded) {
    
    anpi::Matrix<double> a(3,2,4.,anpi::Matrix<double>::Padded);
    ASSERT_EQ(a.rows() , 3) << "Wrong rows";
    ASSERT_EQ(a.cols() , 2) << "Wrong cols";
    ASSERT_EQ(a(0,0) , 4.) << "Wrong value";
    
}

//Initializer list have Padded by default
TEST(MatrixConstructors, PaddedInitializerList) {
    
    anpi::Matrix<float> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    ASSERT_EQ(a.rows() , 3) << "Wrong rows";
    ASSERT_EQ(a.cols() , 4) << "Wrong cols";
    ASSERT_EQ(a(0,0) , 1.f) << "Wrong value";
    ASSERT_EQ(a(1,2) , 7.f) << "Wrong value";
    ASSERT_EQ(a(2,3) , 12.f) << "Wrong value";
    
}

TEST(MatrixConstructors, PaddedCopyConstructor) {
    
    anpi::Matrix<float> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    anpi::Matrix<float> b(a);    
    
    ASSERT_EQ(a, b) << "Not equal";
    ASSERT_EQ(b.rows() , 3) << "Wrong rows";
    ASSERT_EQ(b.cols() , 4) << "Wrong cols";
    ASSERT_NE(b.data() , a.data()) << "Wrong data";
    
}

TEST(MatrixConstructors, PaddedMoveConstructor) {
    
    anpi::Matrix<float> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    anpi::Matrix<float> b(std::move(a));   
    
    ASSERT_EQ(b.rows() , 3) << "Wrong rows";
    ASSERT_EQ(b.cols() , 4) << "Wrong cols";
    ASSERT_TRUE(a.empty()) << "Not moved";
    
}

TEST(MatrixConstructors, PaddedMemConstructor) {
    
    anpi::Matrix<float> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    anpi::Matrix<float> b(a.rows(),a.cols(),a.data());
    
    ASSERT_EQ(a , b) << "Not equal";
    ASSERT_EQ(b.rows() , 3) << "Wrong rows";
    ASSERT_EQ(b.cols() , 4) << "Wrong cols";
    ASSERT_NE(b.data(), a.data()) << "Wrong data";
    
}

//Functions

TEST(MatrixFunctions, MoveAssigment) {
    
    anpi::Matrix<int> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    anpi::Matrix<int> c(a);
    anpi::Matrix<int> b;
    b=std::move(a);
    
    ASSERT_TRUE(a.empty()) << "Not moved";
    ASSERT_FALSE(b.empty()) << "Not moved";
    ASSERT_EQ(b.rows() , 3) << "Wrong rows";
    ASSERT_EQ(b.cols() , 4) << "Wrong cols";
    ASSERT_EQ(b, c) << "Wrong data";
    
}

TEST(MatrixFunctions, Assigment) {
    
    anpi::Matrix<int> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    anpi::Matrix<int> b;
    b=a;
    ASSERT_EQ(a, b) << "Wrong content";
    
}

TEST(MatrixFunctions, Swap) {
    
    anpi::Matrix<int> a = { {1,2,3,4},{5,6,7,8},{9,10,11,12} };
    anpi::Matrix<int> b = { {13,14},{15,16} };

    anpi::Matrix<int> c(a);
    anpi::Matrix<int> d(b);

    ASSERT_EQ(a, c) << "Wrong content";
    ASSERT_EQ(d, b) << "Wrong content";
    c.swap(d);
    ASSERT_EQ(a, d) << "Wrong content";
    ASSERT_EQ(b, c) << "Wrong content";
    
}

//Operators
TEST(OperatorsWithPaddedMatrix, EqualityOperators){
    
    //Initializer list is padded by default
    anpi::Matrix<int> a = {{1,2,3,4},{5,6,7,8},{9,10,11,12}};
    anpi::Matrix<int> b = {{1,2,3,4},{5,6,8,8},{9,10,11,12}};
    
    ASSERT_NE(a , b) << "Equal";
    b(1,2)=7;
    ASSERT_EQ(a , b) << "Not equal"; 
}

TEST(OperatorsWithPaddedMatrix, AddOperator){
    
    //Initializer list is padded by default
    anpi::Matrix<int> a = {{1,2,3},{ 4, 5, 6}};
    anpi::Matrix<int> b = {{7,8,9},{10,11,12}};
    anpi::Matrix<int> r = {{8,10,12},{14,16,18}};
    
    anpi::Matrix<int> c(a);
    c +=  b;
    ASSERT_EQ(c , r) << "Wrong inplace sum"; 
    c=a+b;
    ASSERT_EQ(c , r) << "Wrong external sum"; 


    c=anpi::Matrix<int>{ {1,2,3},{ 4, 5, 6} } + b;
    ASSERT_EQ(c , r) << "Wrong cast sum"; 

    c=a+anpi::Matrix<int>{ {7,8,9},{10,11,12} };
    ASSERT_EQ(c , r) << "Wrong cast sum"; 

}

TEST(OperatorsWithPaddedMatrix, SubstractOperator){
    
    //Initializer list is padded by default
    anpi::Matrix<int> a = { {1,2,3},{ 4, 5, 6} };
    anpi::Matrix<int> b = { {7,8,9},{10,11,12} };
    anpi::Matrix<int> r = { {-6,-6,-6},{-6,-6,-6} };
    
    anpi::Matrix<int> c(a);
    c-=b;
    ASSERT_EQ(c , r) << "Wrong inplace substract"; 
    c=a-b;
    ASSERT_EQ(c , r) << "Wrong external substract"; 


    c=anpi::Matrix<int>{ {1,2,3},{ 4, 5, 6} } - b;
    ASSERT_EQ(c , r) << "Wrong cast substract"; 

    c=a-anpi::Matrix<int>{ {7,8,9},{10,11,12} };
    ASSERT_EQ(c , r) << "Wrong cast substract"; 

}

//+ - Operators with each Type
TEST(OperatorsWithEachType, AddOperatorInt){
    
    //Initializer list is padded by default
    anpi::Matrix<int> a = {{1,2,3},{ 4, 5, 6}};
    anpi::Matrix<int> b = {{7,8,9},{10,11,12}};
    anpi::Matrix<int> r = {{8,10,12},{14,16,18}};
    
    anpi::Matrix<int> c(a);
    c=a+b;
    ASSERT_EQ(c , r) << "Wrong external sum"; 

}

TEST(OperatorsWithEachType, AddOperatorFloat){
    
    //Initializer list is padded by default
    anpi::Matrix<float> a = {{1,2,3},{ 4, 5, 6}};
    anpi::Matrix<float> b = {{7,8,9},{10,11,12}};
    anpi::Matrix<float> r = {{8,10,12},{14,16,18}};
    
    anpi::Matrix<float> c(a);
    c=a+b;
    ASSERT_EQ(c , r) << "Wrong external sum"; 

}

TEST(OperatorsWithEachType, AddOperatorDouble){
    
    //Initializer list is padded by default
    anpi::Matrix<double> a = {{1,2,3},{ 4, 5, 6}};
    anpi::Matrix<double> b = {{7,8,9},{10,11,12}};
    anpi::Matrix<double> r = {{8,10,12},{14,16,18}};
    
    anpi::Matrix<double> c(a);
    c=a+b;
    ASSERT_EQ(c , r) << "Wrong external sum"; 

}

TEST(OperatorsWithEachType, SubstractOperatorInt){
    
    //Initializer list is padded by default
    anpi::Matrix<int> a = {{1,2,3},{ 4, 5, 6}};
    anpi::Matrix<int> b = {{7,8,9},{10,11,12}};
    anpi::Matrix<int> r = { {-6,-6,-6},{-6,-6,-6} };
    
    anpi::Matrix<int> c(a);
    c=a-b;
    ASSERT_EQ(c , r) << "Wrong external sum"; 

}

TEST(OperatorsWithEachType, SubstractOperatorFloat){
    
    //Initializer list is padded by default
    anpi::Matrix<float> a = {{1,2,3},{ 4, 5, 6}};
    anpi::Matrix<float> b = {{7,8,9},{10,11,12}};
    anpi::Matrix<float> r = { {-6,-6,-6},{-6,-6,-6} };
    
    anpi::Matrix<float> c(a);
    c=a-b;
    ASSERT_EQ(c , r) << "Wrong external sum"; 

}

TEST(OperatorsWithEachType, SubstractOperatorDouble){
    
    //Initializer list is padded by default
    anpi::Matrix<double> a = {{1,2,3},{ 4, 5, 6}};
    anpi::Matrix<double> b = {{7,8,9},{10,11,12}};
    anpi::Matrix<double> r = { {-6,-6,-6},{-6,-6,-6} };
    
    anpi::Matrix<double> c(a);
    c=a-b;
    ASSERT_EQ(c , r) << "Wrong external sum"; 

}

//Index Map Tests

TEST(IndexMap, NodesToCurrentIndex){
    
    IndexMap * indexMap = new IndexMap(4,4);
    //First current
    int i0 = indexMap->getXFromNodes(0,0,0,1);//0,0 to 0,1
    //Random current
    int i19 = indexMap->getXFromNodes(2,2,3,2);//2,2 to 3,2
    //Last current
    int i23 = indexMap->getXFromNodes(3,2,3,3);// 3,2 to 3,3

    ASSERT_EQ(i0 , 0) << "Wrong nodes to current"; 
    ASSERT_EQ(i19 , 19) << "Wrong nodes to current"; 
    ASSERT_EQ(i23 , 23) << "Wrong nodes to current"; 

}

TEST(IndexMap, CurrentIndexToNode){
    
    IndexMap * indexMap = new IndexMap(4,4);
    
    //Random current
    int r1 = 0;
    int c1 = 0;
    int r2 = 0;
    int c2 = 0;

    indexMap->getNodesFromX(r1,c1,r2,c2,19); //Nodes of current 19
   

    ASSERT_EQ(r1 , 2) << "Wrong nodes to current"; 
    ASSERT_EQ(c1 , 2) << "Wrong nodes to current"; 
    ASSERT_EQ(r2 , 3) << "Wrong nodes to current"; 
    ASSERT_EQ(c2 , 2) << "Wrong nodes to current"; 

}

TEST(IndexMap, TwoSidesIntegrity){
    
    IndexMap * indexMap = new IndexMap(4,4);
    
    //Random current
    int r1 = 0;
    int c1 = 0;
    int r2 = 0;
    int c2 = 0;

    indexMap->getNodesFromX(r1,c1,r2,c2,3); //Nodes of current 3
    int i = indexMap->getXFromNodes(r1,c1,r2,c2);
   

    ASSERT_EQ(i , 3) << "Wrong nodes to current"; 


}
