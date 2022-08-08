// CMSC 341 - Spring 2022 - Project 4
//dnadb.cpp
//gurnoor Hira
#include "dnadb.h"
DnaDb::DnaDb(int size, hash_fn hash){
    if(size < MINPRIME){
        size = MINPRIME;
    }
    if(size > MAXPRIME){
        size = MAXPRIME;
    }
    if(!isPrime(size)){
        size = findNextPrime(size);
    }
    m_currentlyHashing = false;
    m_oldNumDeleted = 0;
    m_currNumDeleted = 0;
    m_oldSize = 0;
    m_currentSize = 0;
    m_oldCap = 0;
    m_currentCap = size;
    m_oldTable = nullptr;
    m_currentTable = new DNA[m_currentCap];
    m_hash = hash;
}

DnaDb::~DnaDb(){
    m_hash = nullptr;
    m_oldSize = 0;
    m_currentSize = 0;
    m_oldCap = 0;
    m_currentCap = 0;
    m_oldNumDeleted = 0;
    m_currNumDeleted = 0;
    m_currentlyHashing = false;
    if(m_oldTable != nullptr){
        delete [] m_oldTable;
    }
    if( m_currentTable != nullptr){
        delete [] m_currentTable;
    }
}

bool DnaDb::insert(DNA dna){
    //inserting into new table
    bool insertComplete = false;
    int quadCounter = 0;
    unsigned int hashNum = m_hash(dna.getSequence());
    hashNum = hashNum % m_currentCap;
    const unsigned int tempNum= hashNum; //holds old hash value
    if(m_currentTable[hashNum].getSequence().empty()){
        m_currentTable[hashNum] = dna;
        m_currentSize++;
        insertComplete = true;
    }
    else{
        //search for open index while quadratic probing
        while(!m_currentTable[hashNum].getSequence().empty() && !(m_currentTable[hashNum] == DELETED)){
            hashNum = tempNum;
            hashNum = (hashNum + (quadCounter * quadCounter)) % m_currentCap;
            quadCounter++;
        }
        m_currentTable[hashNum] = dna;
        m_currentSize++;
        insertComplete = true;
    }
    if(m_currentlyHashing == true){
        reHashOldTable();
    }
    else{
        if(lambda() > .5){
            createNewHashCurrTable();
            m_currentlyHashing = true;
        }
    }
    

    return insertComplete;

}

bool DnaDb::remove(DNA dna){
    bool removed = false;
    if(m_oldTable != nullptr){ //check if old table exists
        for(unsigned int i =0; i<m_oldCap;i++){
            if(!m_oldTable[i].getSequence().empty()){ //if index empty then search for dna
                if(m_oldTable[i] == dna){
                    m_oldTable[i] = DELETED;
                    m_oldNumDeleted++;
                    removed = true;
                }
            }
        }
        if(m_currentlyHashing == true){ //working w old table rehash new table
            reHashCurrTable();
        }
        else if(!m_currentlyHashing){
            if(oldDeletedRatio() > .8){
                m_currentlyHashing = true;
                createNewHashOldTable();
            }
        }
    }
    
    if(m_currentTable != nullptr){ //check if current table exist and do same process
        for(unsigned int i = 0;i<m_currentCap;i++){
            if(!m_currentTable[i].getSequence().empty()){
                if(m_currentTable[i] == dna){
                    m_currentTable[i] = DELETED;
                    m_currNumDeleted++;
                    removed = true;
                }
            }
        }
        if(m_currentlyHashing == true){
            reHashOldTable();
        }
        else if(!m_currentlyHashing){
            if(deletedRatio() > .8){
                m_currentlyHashing = true;
                createNewHashCurrTable();
            }
        }
        
    }
    return removed;
}

DNA DnaDb::getDNA(string sequence, int location){
    DNA nullDna; //null dna to return if dna not found
    if(m_oldTable != nullptr){ //search old table 
        for(unsigned int i = 0;i < m_oldCap; i++){
            if(!m_oldTable[i].getSequence().empty()){ //if index not empty then search it
                if(m_oldTable[i].getSequence() == sequence && m_oldTable[i].getLocId() == location){
                    return m_oldTable[i];
                }
            }
        }
    }
    if(m_currentTable != nullptr){ //repeat process for current table
        for(unsigned int i; i <m_currentCap; i++){
            if(!m_currentTable[i].getSequence().empty()){
                if(m_currentTable[i].getSequence() == sequence && m_currentTable[i].getLocId() == location){
                    return m_currentTable[i];
                }
            }
        }
    }
    cout << "DNA not found" << endl;
    return nullDna;
    

}

float DnaDb::lambda() const {
    float size = float(m_currentSize);
    float capacity = float(m_currentCap);
    return size/capacity; //returns float form of lambda
      
}

float DnaDb::deletedRatio() const { //deleted ratio of current table
    float numDelete = 0.0;
    for(unsigned int i = 0; i <m_currentCap; i++){ //count deleted keys
        if(!m_currentTable[i].getSequence().empty()){
            if(m_currentTable[i].getSequence() == DELETEDKEY){
                numDelete += 1.0;
            }
        }
    }
    float capacity = float(m_currentSize); //deleted key /occupied buckets
    return numDelete / capacity;
}

void DnaDb::dump() const {
    cout << "Dump for current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

bool DnaDb::isPrime(int number){
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int DnaDb::findNextPrime(int current){
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) { 
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0) 
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}

DNA::DNA(string sequence, int location) {
    if ((location >= MINLOCID && location <= MAXLOCID) ||
        (location == 0 && sequence == "DELETED")){
        // this is a normal or a DELETED object
        m_sequence = sequence;
        m_location = location;
    }
    else{
        // this is the empty object
        m_sequence = "";
        m_location = 0;
    }
}

string DNA::getSequence() const {
    return m_sequence;
}

int DNA::getLocId() const {
    return m_location;
}

// Overloaded assignment operator
const DNA& DNA::operator=(const DNA& rhs){
    if (this != &rhs){
        m_sequence = rhs.m_sequence;
        m_location = rhs.m_location;
    }
    return *this;
}

// Overloaded insertion operator.  Prints DNA's sequence (key),
// and the location ID. This is a friend function in DNA class.
ostream& operator<<(ostream& sout, const DNA &dna ) {
    if (!dna.m_sequence.empty())
        sout << dna.m_sequence << " (Location ID " << dna.m_location << ")";
    else
        sout << "";
  return sout;
}

// Overloaded equality operator. This is a friend function in DNA class.
// To test inequality we may negate the results of this operator.
bool operator==(const DNA& lhs, const DNA& rhs){
    return ((lhs.m_sequence == rhs.m_sequence) && (lhs.m_location == rhs.m_location));
}
void DnaDb::reHashOldTable(){
    int numD = m_oldSize * .25; //25% of  size
    unsigned int i = 0;
    int currD = 0;
    //keep going till 25% is hit or iterate through full table
    while(i<m_oldCap && currD <= numD){
        if(!m_oldTable[i].getSequence().empty()){
            if(m_oldTable[i].getSequence() != DELETEDKEY){
                //temp dna to deep copy on new table
                DNA temp(m_oldTable[i].getSequence(),m_oldTable[i].getLocId());
                insertHelperOldTable(temp);
                m_oldTable[i] = DELETED;
                currD++;
            }
        }
        i++;
    }
    //rehash complete
    if(i == m_oldCap){
        m_currentlyHashing = false;
        delete [] m_oldTable;
        m_oldSize = 0;
        m_oldNumDeleted = 0;
        m_oldTable = nullptr;
        m_oldCap = 0;
    }
}
void DnaDb::createNewHashCurrTable(){
    if(m_oldTable != nullptr){ //if not cleared then clear
        delete [] m_oldTable;
        m_oldTable = nullptr;
    }
    m_oldCap = ((m_currentSize = m_currNumDeleted) * 4);
    //find prime num for new cap
    if(!isPrime(m_oldCap)){
        m_oldCap = findNextPrime(m_oldCap);
    }
    //start hash and allocate memory
    cout <<"rehashing" << endl;
    m_oldTable = new DNA[m_oldCap];
    m_currentlyHashing = true;
         

}
void DnaDb::insertHelperOldTable(DNA dna){
    int m_quadCounter = 0;
    unsigned int hashNum = m_hash(dna.getSequence());
    hashNum = hashNum % m_oldCap;
    if(m_oldTable[hashNum].getSequence().empty()){
        m_oldTable[hashNum] = dna;
        m_oldSize++;
    }
    else{
        while(!m_oldTable[hashNum].getSequence().empty() && !(m_oldTable[hashNum] == DELETED)){
            hashNum = (hashNum + (m_quadCounter * m_quadCounter)) % m_oldCap;
            m_quadCounter++;
        }
        m_oldTable[hashNum] = dna;
        m_oldSize++;
    }

}
void DnaDb::insertHelperCurrTable(DNA dna){
    int m_quadCounter = 0;
    unsigned int hashNum = m_hash(dna.getSequence());
    hashNum = hashNum % m_currentCap;
    if(m_currentTable[hashNum].getSequence().empty()){
        hashNum = (hashNum + (m_quadCounter * m_quadCounter)) % m_currentCap;
        m_currentSize++;
    }
    else{
        while(!m_currentTable[hashNum].getSequence().empty() && !(m_currentTable[hashNum] == DELETED)){
            hashNum = (hashNum + (m_quadCounter * m_quadCounter)) % m_currentCap;
            m_quadCounter++;
        }
        m_currentTable[hashNum] = dna;
        m_currentSize++;
    }

}

float DnaDb::oldDeletedRatio() const{
    float numDelete = 0.0;
    for(unsigned int i = 0; i <m_oldCap; i++){
        if(!m_oldTable[i].getSequence().empty()){
            if(m_oldTable[i].getSequence() == DELETEDKEY){
                numDelete += 1.0;
            }
        }
    }
    float capacity = float(m_oldSize);
    return numDelete / capacity;
}
void DnaDb::reHashCurrTable(){
    int numD = m_currentSize * .25;
    unsigned int i = 0;
    int currD = 0;
    while( i < m_currentCap && currD <= numD){
        if(!m_currentTable[i].getSequence().empty()){
            if(m_currentTable[i].getSequence() != DELETEDKEY){
                DNA temp(m_currentTable[i].getSequence(),m_currentTable[i].getLocId());
                insertHelperCurrTable(temp);
                m_currentTable[i] = DELETED;
                currD++;
            }
        }
        i++;
    }
    //rehash complete
    if(i == m_currentCap){
        m_currNumDeleted = 0;
        m_currentlyHashing = false;
        m_currentSize = 0;
        delete [] m_currentTable;
        m_currentTable = nullptr;
        m_currentCap = 0;
    }
}
void DnaDb::createNewHashOldTable(){
        if(m_currentTable != nullptr){ //if not cleared then clear
            delete [] m_currentTable;
            m_currentTable = nullptr;
        }
        m_currentCap = ((m_oldSize = m_oldNumDeleted) * 4);
        if(!isPrime(m_currentCap)){
            m_currentCap = findNextPrime(m_currentCap);
        }
        cout <<"rehashing" << endl;
        m_currentTable = new DNA[m_oldCap];
        m_currentlyHashing = true;   
}