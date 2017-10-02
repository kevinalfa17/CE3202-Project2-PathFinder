/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 */



namespace anpi
{
  // ------------------------
  // Implementation of Matrix
  // ------------------------

  template<typename T>
  Matrix<T>::Matrix() : _data(0),_rows(0),_cols(0),_dcols(0) {}
  
//Initialization with type
#if defined(IS_SIMD_ACTIVE) && defined(IS_SIMD_AVAILABLE) 
  //DoNotInitialize and Padded possible
  template<typename T>
  Matrix<T>::Matrix(const size_t r, const size_t c, const InitializationType): _data(0),_rows(0),_cols(0),_dcols(0) {
	  switch(InitializationType){
	  case DoNotInitialize:
		  if(() 
		  if(((sizeof(T)*c)/16) != 1 || (r*c)%4 != 0){
			  throw PaddingWarningException();
		  }
		  
		  allocate(r,c);
		  break;
	  case Padded:
		  padded_allocate(r,c);
		  break;
	  }
	  
  }
#else
  //Only DoNotInitialize is possible
  template<typename T>
  Matrix<T>::Matrix(const size_t r, const size_t c, const InitializationType): _data(0),_rows(0),_cols(0),_dcols(0) {
    allocate(r,c);
  }
#endif
  
  
  //Initialization with initial value
  template<typename T>
  Matrix<T>::Matrix(const size_t r, const size_t c, const T initVal):Matrix(r,c,DoNotInitialize) {
    fill(initVal);
  }
  
  //Initialization with  given pointer
  template<typename T>
  Matrix<T>::Matrix(const size_t r,  const size_t c, const T *const initMem): Matrix(r,c,DoNotInitialize) {
    fill(initMem);
  }

  //Initialization with initializer list
  template<typename T>
  Matrix<T>::Matrix(std::initializer_list< std::initializer_list<T> > lst):Matrix(lst.size(),(lst.size()>0) ? lst.begin()->size() : 0,DoNotInitialize)
  {
    size_t j = 0;
    for (const auto& r : lst) {
      // Verify that all rows have the same number of columns
      assert( (r.size()==_cols) && "Check number of columns in each row");
      size_t rowIdx=j*_cols;
      for (const auto& c : r) {
	_data[rowIdx++] = c;
      }
      ++j;
    }
  }
    
  //Copy constructor (Deep copy)
  template<typename T>
  Matrix<T>::Matrix(const Matrix<T>& other): Matrix(other._rows,other._cols,DoNotInitialize) {
    fill(other.data());
  }
  
  //Transfer constructor 
  template<typename T>
  Matrix<T>::Matrix(Matrix<T>&& other): _data(other._data),_rows(other._rows),_cols(other._cols), _dcols(other._dcols) {
    other._data=0; // remove references to the data in other
    other._rows=0; 
    other._cols=0;
    other._dcols=0;
  }
  
  //Relase all memory
  template<typename T>
  Matrix<T>::~Matrix() {
	  deallocate();
	  _data=0;
	  _rows=0;
	  _cols=0;
	  _dcols=0;
  }
  
//Functions with InitilizationType Parameter
#if defined(IS_SIMD_ACTIVE) && defined(IS_SIMD_AVAILABLE) 
  
  //Initialization with initial value and initialization type
  template<typename T>
  Matrix<T>::Matrix(const size_t r, const size_t c, const T initVal, const InitializationType):Matrix(r,c,InitializationType) {
	  
	  switch(InitializationType){
	  case DoNotInitialize:
		  fill(initVal);
		  break;
	  case Padded:
		  padded_fill(initVal);
		  break;
	  }  
  }

  //Initialization with  given pointer
  template<typename T>
  Matrix<T>::Matrix(const size_t r,  const size_t c, const T *const initMem,const InitializationType): Matrix(r,c,InitializationType) {
	  switch(InitializationType){
	  case DoNotInitialize:
		  fill(initMem);
		  break;
	  case Padded:
		  padded_fill(initMem);
		  break;
	  }
  }
 
  //Copy constructor (Deep copy)
  template<typename T>
  Matrix<T>::Matrix(const Matrix<T>& other, const InitializationType): Matrix(other._rows,other._cols,InitializationType) {
	  
	  switch(InitializationType){
	  case DoNotInitialize:
		  fill(other.data());
		  break;
	  case Padded:
		  padded_fill(other.data());
		  break;
	  }
  }
  
#endif

 

  template<typename T>
  void Matrix<T>::swap(Matrix<T>& other) {
    std::swap(_data,other._data);
    std::swap(_rows,other._rows);
    std::swap(_cols,other._cols);
  }
  
#if defined(IS_SIMD_ACTIVE) && defined(IS_SIMD_AVAILABLE) //SIMD ALLOCATION
  template<typename T>
  void Matrix<T>::allocate(const size_t r, const size_t c) {
	  
	  // Only reserve if the desired size is different to the current one
	  if ( (r!=_rows) || (c!=_cols) ) {
		  deallocate(); //Free aligned memory
		  const size_t total = r*c; // Total number of required entries
		  T* alignData = (T*) _mm_malloc (total * sizeof(T), 16);//Align total T spaces to 16-byte for SSE

		  _data = (total > 0) ? alignData : 0;

		  _rows=r;
		  _cols=c;
		  _dcols=c;
	  }
  }
  
  void Matrix<T>::padded_allocate(const size_t r, const size_t c) {
	  
	  // Only reserve if the desired size is different to the current one
	  if ( (r!=_rows) || (c!=_cols) ) {
		  deallocate(); //Free aligned memory
		  
		  int res = (sizeof(T)*c)/16; //Aux variable to check padding required 
		  int cols_per_row = 16/sizeof(T); //Max number of cols for T type
				  
		  if(res > 1){ //Extra column padding required
			  _dcols = (c\cols_per_row)*cols_per_row+cols_per_row; //Pay attention that \ instead od / is used!		  
		  }
		  else if(res < 1){ //Minimum colums 16/sizeof(T)required. Ex: 4 in float or 2 in double
			  _dcols = cols_per_row;
		  }
		  else{ //No padding required (cols are  multiple of 16)
			  _dcols = c;  
		  }
		  
		  const size_t total = r*_dcols; // Total number of required entries

		  
		  T* alignData = (T*) _mm_malloc (total * sizeof(T), 16);//Align total*mult T spaces to 16-byte for SSE

		  _data = (total > 0) ? alignData : 0;
		  _rows=r;
		  _cols=c;

	  }
  }
#else //Standar allocation
  template<typename T>
  void Matrix<T>::allocate(const size_t r, const size_t c) {
	  // only reserve if the desired size is different to the current one
	  if ( (r!=_rows) || (c!=_cols) ) {
		  deallocate(); //Free aligned memory
		  const size_t total = r*c; // total number of required entries
		  _data = (total > 0) ? new T[total] : 0;

		  _rows=r;
		  _cols=c;
		  _dcols=c;
	  }
  }
#endif

//Fill function
#if defined(IS_SIMD_ACTIVE) && defined(IS_SIMD_AVAILABLE) //SIMD fill
  //Float fill
  template<>
  void Matrix<float>::fill(const float val) {

	  int SSELength = (_rows*_cols)/4; //fill 4 values at the same time
	  __m128 * dataSSE =  (__m128*) _data;

	  for(int i = 0; i < SSELength; i++){
		  dataSSE[i] = _mm_set1_ps(val); // dataSSE[i] = [val,val,val,val]
	  }

  }

  //Double fill
  template<>
  void Matrix<double>::fill(const double val) {

	  int SSELength = (_rows*_cols)/2; //fill 2 values at the same time
	  __m128d * dataSSE =  (__m128d*) _data;

	  for(int i = 0; i < SSELength; i++){
		  dataSSE[i] = _mm_set1_pd(val); // dataSSE[i] = [val,val]
	  }

  }

  //int fill
  template<>
  void Matrix<int>::fill(const int val) {

	  int SSELength = (_rows*_cols)/2; //fill 4 values at the same time
	  __m128i * dataSSE =  (__m128i*) _data;

	  for(int i = 0; i < SSELength; i++){
		  dataSSE[i] = _mm_set1_epi32(val); // dataSSE[i] = [val,val,val,val]
	  }

  }

  //Type is different to float, int or double, then do the standard fill
  template<typename T>
   void Matrix<T>::fill(const T val) {
 	  T* end = _data + (_rows*_cols);
 	  for (T* ptr = _data;ptr!=end;++ptr) {
 		  *ptr = val;
 	  }
   }
  
  
  //Float padded_fill
   template<>
   void Matrix<float>::padded_fill(const float val) {

 	  int SSELength = (_rows*_dcols)/4; //fill 4 values at the same time
 	  __m128 * dataSSE =  (__m128*) _data;

 	  for(int i = 0; i < SSELength; i++){
 		  dataSSE[i] = _mm_set1_ps(val); // dataSSE[i] = [val,val,val,val]
 	  }

   }

   //Double padded_fill
   template<>
   void Matrix<double>::padded_fill(const double val) {

 	  int SSELength = (_rows*_dcols)/2; //fill 2 values at the same time
 	  __m128d * dataSSE =  (__m128d*) _data;

 	  for(int i = 0; i < SSELength; i++){
 		  dataSSE[i] = _mm_set1_pd(val); // dataSSE[i] = [val,val]
 	  }

   }

   //int padded_fill
   template<>
   void Matrix<int>::padded_fill(const int val) {

 	  int SSELength = (_rows*_dcols)/2; //fill 4 values at the same time
 	  __m128i * dataSSE =  (__m128i*) _data;

 	  for(int i = 0; i < SSELength; i++){
 		  dataSSE[i] = _mm_set1_epi32(val); // dataSSE[i] = [val,val,val,val]
 	  }

   }

   //Type is different to float, int or double, then do the standard fill
   template<typename T>
    void Matrix<T>::padded_fill(const T val) {
  	  T* end = _data + (_rows*_dcols);
  	  for (T* ptr = _data;ptr!=end;++ptr) {
  		  *ptr = val;
  	  }
    }
   
   //Memcpy padded
   template<typename T>
   void Matrix<T>::padded_fill(const T* mem) {
     std::memcpy(_data,mem,sizeof(T)*_rows*_dcols);
   }
   
#else //Standard fill
  
  template<typename T>
  void Matrix<T>::fill(const T val) {
	  T* end = _data + (_rows*_cols);
	  for (T* ptr = _data;ptr!=end;++ptr) {
		  *ptr = val;
	  }
  }
  

#endif
  
  template<typename T>
  void Matrix<T>::fill(const T* mem) {
    std::memcpy(_data,mem,sizeof(T)*_rows*_cols);
  }
  
  
  //Assign and comparisson operators
  template<typename T>
   Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) {
 	  allocate(other._rows,other._cols);
 	  fill(other.data());
 	  return *this;
   }

   // = Operator
   template<typename T>
   Matrix<T>& Matrix<T>::operator=(Matrix<T>&& other) {
 	  if (_data!=other._data) { // alias detection first
 		  deallocate();

 		  _data=other._data;
 		  _rows=other._rows;
 		  _cols=other._cols;
 		  _dcols=other._dcols;
 		  

 		  other._data=0;
 		  other._rows=0;
 		  other._cols=0;
 		  other._dcols=0;
 	  }

 	  return *this;
   }
   
   template<typename T>
   bool Matrix<T>::operator==(const Matrix<T>& other) const {
     if (&other==this) return true; // alias detection

     // same size of matrices?
     if ((other._rows != _rows) ||
 	(other._cols != _cols)) return false;

     // check the content with pointers
     return (memcmp(_data,other._data,entries()*sizeof(T))==0);
   }

   template<typename T>
   bool Matrix<T>::operator!=(const Matrix<T>& other) const {
     if (&other==this) return false; // alias detection

     // same size of matrices?
     if ((other._rows != _rows) ||
 	(other._cols != _cols)) return true;

     // check the content with pointers
     return (memcmp(_data,other._data,entries()*sizeof(T))!=0);
   }

  // ARITHMETIC OPERATORS
  template<typename T>
  Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& other) {

    assert( (other._rows==_rows) && (other._cols==_cols) );

    T* here        = _data;
    T *const end   = _data + entries();
    const T* there = other._data;

    for (;here!=end;) {
      *here++ += *there++;
    }

    return *this;
  }

  template<typename T>
  Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& other) {

    assert( (other._rows==_rows) && (other._cols==_cols) );

    T* here        = _data;
    T *const end   = _data + entries();
    const T* there = other._data;

    for (;here!=end;) {
      *here++ -= *there++;
    }

    return *this;
  }
  

  
  template<class T>
  Matrix<T> operator+(const Matrix<T>& a,
		      const Matrix<T>& b) {
    
    assert( (a.rows()==b.rows()) && (a.cols()==b.cols()) );
    
    Matrix<T> c(a.rows(),a.cols(),Matrix<T>::DoNotInitialize);
    T* cptr        = c.data();
    T *const end   = cptr + a.entries();

    const T* aptr = a.data();
    const T* bptr = b.data();

    for (;cptr!=end;) {
      *cptr++ = *aptr++ + *bptr++;
    }

    return c;
  }

  template<class T>
  Matrix<T> operator-(const Matrix<T>& a,
		      const Matrix<T>& b) {
    
    assert( (a.rows()==b.rows()) && (a.cols()==b.cols()) );
    
    Matrix<T> c(a.rows(),a.cols(),Matrix<T>::DoNotInitialize);
    T* cptr        = c.data();
    T *const end   = cptr + a.entries();

    const T* aptr = a.data();
    const T* bptr = b.data();

    for (;cptr!=end;) {
      *cptr++ = *aptr++ - *bptr++;
    }

    return c;
  }
  
  
  template<class T>
  Matrix<T> operator*(const Matrix<T>& a, const Matrix<T>& b) {

	  Matrix<T> c(a.rows(),b.cols(),T(0));
	  
	  if(a.cols() != b.rows()){
		  throw MatrixException();
		  
	  }
	  else{

		  // Multiplying matrix a and b and storing in array mult.
		  for(int i = 0; i < a.rows(); i++){
			  for(int j = 0; j < b.cols(); j++){
				  for(int k = 0; k < a.rows(); k++)
				  {
					  c[i][j] += a[i][k] * b[k][j];
				  }	
			  }
		  }

	  }

	  return c;

  }

#if defined(IS_SIMD_ACTIVE) && defined(IS_SIMD_AVAILABLE) 
  template<typename T>
  void Matrix<T>::deallocate() {
	  _mm_free ((void *) _data);
  }	
#else //Standard memory liberation
  template<typename T>
  void Matrix<T>::deallocate() {
	  delete[] _data;
  }	
#endif
  
  
} // namespace ANPI
