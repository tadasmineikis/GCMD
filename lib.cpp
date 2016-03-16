#include "lib.h"
#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

using namespace std;


#define   NULL_PTR   (void *) 0
#define   ZERO       1.0E-10
//#define   INFINITY   1.0E15
#define   UL         unsigned long

float LogInterp(float y2, float y1, float x2, float x1, float x)
{
	return pow(10, log10(y1) + (log10(y2) - log10(y1))/(log10(x2) - log10(x1))*( log10(x)-log10(x1) ) );
}

float LinInterp(float y2, float y1, float x2, float x1, float x)
{
	return y1 + (y2 - y1)/(x2 - x1)*(x - x1);
	//return ((x - x1)*y2 + (x2 - x)*y1)/(x2 - x1);
}

double LinInterpP(double y2, double y1, double x2, double x1, double x)
{
	return y1 + (y2 - y1)/(x2 - x1)*(x - x1);
	//return ((x - x1)*y2 + (x2 - x)*y1)/(x2 - x1);
}


  /*
   * The function
   *      void  **matrix()
   * reserves dynamic memory for a two-dimensional matrix
   * using the C++ command new . No initialization of the elements.
   * Input data:
   *  int row      - number of  rows
   *  int col      - number of columns
   *  int num_bytes- number of bytes for each
   *                 element
   * Returns a void  **pointer to the reserved memory location.
   */

void **matrix(int row, int col, int num_bytes)
  {
  int      i, num;
  char     **pointer, *ptr;

  pointer = new(nothrow) char* [row];
  if(!pointer) {
    cout << "Exception handling: Memory allocation failed";
    cout << " for "<< row << "row addresses !" << endl;
    return NULL;
  }
  i = (row * col * num_bytes)/sizeof(char);
  pointer[0] = new(nothrow) char [i];
  memset(pointer[0], 0, i*sizeof(char));
  if(!pointer[0]) {
    cout << "Exception handling: Memory allocation failed";
    cout << " for address to " << i << " characters !" << endl;
    return NULL;
  }
  ptr = pointer[0];
  num = col * num_bytes;
  for(i = 0; i < row; i++, ptr += num )   {
    pointer[i] = ptr;
  }

  return  (void **)pointer;

  } // end: function void **matrix()

    /*
     * The function
     *      void free_matrix()
     * releases the memory reserved by the function matrix()
     *for the two-dimensional matrix[][]
     * Input data:
     *  void far **matr - pointer to the matrix
     */

void free_matrix(void **matr)
{

  delete [] (char *) matr[0];
  delete [] matr;

}  // End:  function free_matrix()

      /*
      ** The function
      **              rk4()
      ** takes a set of variables y[1:n] for the function y(x) together with the
      ** derivatives dydx[1:n] and uses the fourth-order Runge-Kutta method to
      ** advance the solution over an interval h and return incremented variables
      ** as yout[1:n], which not need to be a disstinct arra from y[1:n]. The
      ** users supply the routine derivs(x,y,dydx), which returns the derivatives
      ** dydx at x.
      */
/*int round(double number)
{
    return number < 0.0 ? (int)ceil(number - 0.5) : (int)floor(number + 0.5);
}
*/
