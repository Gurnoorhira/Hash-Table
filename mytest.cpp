// CMSC 341 - Spring 2022 - Project 4
#include "dnadb.h"
#include "dnadb.cpp"
#include <random>
#include <vector>
enum RANDOM {UNIFORMINT, UNIFORMREAL, NORMAL};
class Random {
public:
    Random(int min, int max, RANDOM type=UNIFORMINT, int mean=50, int stdev=20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL){
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            //the mean and standard deviation can change by passing new values to constructor 
            m_normdist = std::normal_distribution<>(mean,stdev);
        }
        else if (type == UNIFORMINT) {
            //the case of UNIFORMINT to generate integer numbers
            // Using a fixed seed value generates always the same sequence
            // of pseudorandom numbers, e.g. reproducing scientific experiments
            // here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min,max);
        }
        else{ //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
        }
    }
    void setSeed(int seedNum){
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }

    int getRandNum(){
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if(m_type == NORMAL){
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while(result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT){
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum(){
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result*100.0)/100.0;
        return result;
    }
    
    private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;//normal distribution
    std::uniform_int_distribution<> m_unidist;//integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal;//real uniform distribution

};
class Tester{
    public:
    void testNonColliding();
    void testColliding();
    private:
};

unsigned int hashCode(const string str);
string sequencer(int size, int seedNum);

int main(){
    Tester test;
    test.testNonColliding();
    test.testColliding();

    return 0;
}
unsigned int hashCode(const string str) {
   unsigned int val = 0 ;
   const unsigned int thirtyThree = 33 ;  // magic number from textbook
   for ( int i = 0 ; i < str.length(); i++)
      val = val * thirtyThree + str[i] ;
   return val ;
}
string sequencer(int size, int seedNum){
    //this function returns a random DNA sequence
    string sequence = "";
    Random rndObject(0,3);
    rndObject.setSeed(seedNum);
    for (int i=0;i<size;i++){
        sequence = sequence + ALPHA[rndObject.getRandNum()];
    }
    return sequence;
}
void Tester::testNonColliding(){
    Random RndLocation(MINLOCID,MAXLOCID);
    int tempLoc[20] = {0};
    DnaDb dnadb(MINPRIME, hashCode);
    int temp = 0;
    for(int i = 0; i < 20; i++){
        temp = RndLocation.getRandNum();
        tempLoc[i] = temp;
    }

    cout << endl;
    cout <<"test insert non collision" << endl;
    DNA testFirst("first insertion",tempLoc[1]);
    DNA testSecond("second insertion",tempLoc[3]);
    DNA testThird("third insertion",tempLoc[6]);
    dnadb.insert(testFirst);
    dnadb.insert(DNA("fourth test",tempLoc[7]));
    dnadb.insert(DNA("fifth test",tempLoc[8]));
    dnadb.insert(testSecond);
    dnadb.insert(testThird);
    dnadb.insert(DNA("sixth test",tempLoc[2]));
    dnadb.insert(DNA("seventh test",tempLoc[4]));
    unsigned int hashNum = hashCode("first insertion");
    hashNum = hashNum % MINPRIME;
    cout << "first insertion should be inserted at " << hashNum << endl;
    cout << endl;

    hashNum = hashCode("fourth test");
    hashNum = hashNum % MINPRIME;
    cout << "fourth test should be inserted at " << hashNum << endl;
    cout << endl;
    cout << "dumping insertion of non collision" << endl;
    dnadb.dump();
    cout << endl;

    cout <<"check if size updating" << endl;
    cout <<"num items = " << dnadb.m_currentSize << endl;
    cout << endl;

    cout << "Test get, non collision" << endl;
    DNA dna2;
    dna2 = dnadb.getDNA("first insertion",tempLoc[1]);
    if(!(dna2.getSequence() == "" || dna2.getLocId() == 1)){
        cout << "first insertion found" << endl;
        cout << endl;
    }

    cout << "test removal non collision" << endl;
    cout << "removing " << testSecond.getSequence() << " ID = " << testSecond.getLocId() << endl;
    if(dnadb.remove(testSecond)){
        cout <<"remove complete" << endl;
    }
    cout << endl;
    cout << "removing " << testThird.getSequence() << " ID = " << testThird.getLocId() << endl;
    if(dnadb.remove(testThird)){
        cout <<"remove complete" << endl;
    }
    cout <<"dump after removals" << endl;
    dnadb.dump();
}
void Tester::testColliding(){
    cout << "Testing insertion, collision without rehash" << endl;
    cout << endl;
    Random RndLocation(MINLOCID,MAXLOCID);
    int tempLoc[50] = {0};
    DnaDb dnadb(MINPRIME, hashCode);
    int temp = 0;
    int sec = 0;
    for(int i = 0; i < 50; i++){
        temp = RndLocation.getRandNum();
        if(i % 3 == 0){
            tempLoc[sec] = temp;
            sec++;
        }
        if(i %3 != 0){
            dnadb.insert(DNA("first insertion",temp));
        }
        else{
            dnadb.insert(DNA("second insertion",temp));
        }
    }
    cout<< "dump after 50 insertions with minprime buckets" <<endl;
    dnadb.dump();
    cout<<endl;
    cout << "Test removal collisiding keys without trigerring rehash" << endl;
    for(int i = 0; i<14;i++){
        dnadb.remove(DNA("second insertion",tempLoc[i]));
    }
    cout << "dump after removing 14 " << endl;
    dnadb.dump();
    
    cout << endl;
    cout <<"Test insert of colliding with rehash" << endl;
    Random RndLocation2(MINLOCID,MAXLOCID);
    int tempLoc2[100] = {0};
    DnaDb rehashingTable(MINPRIME, hashCode);
    int temp2 = 0;
    int sec2 = 0;
    cout << "inserting 100 buckets" << endl;
    for(int i = 0; i <100; i++){
        temp2 = RndLocation2.getRandNum();
        tempLoc2[sec2] = temp2;
        sec2++;
        rehashingTable.insert(DNA("first insertion",temp2));
    }
    cout <<"dump after insertion and rehash" << endl;
    rehashingTable.dump();
    cout << endl;
    cout <<"test removal of colliding with rehash" << endl;
    for(int i = 0;i < 84; i++){
        rehashingTable.remove(DNA("first insertion",tempLoc[i]));
    }
    cout << endl;
    cout << "dump after remove and rehash" << endl;
    rehashingTable.dump();

}
