#include <fstream>
#include <cmath>
//#include <limits.h>

#include "util.h"



//////////
// indexx
void indexx::operator() (const int n, const vector<double>& arrin, 
			 vector<int>& indx)
{
	// Returns an indexed array according to arrin
	int i, j, L, ir, indxt;
	double q;	//q and arrin are doubles!!!
	if(n<=0)
		return;
	for(j=0; j<n; j++) indx[j] = j;
	L = n/2;
	ir = n-1;
	for(;;) {
		if(L>0) {
			indxt = indx[--L];
			q = arrin[indxt];
		}
		else {
			indxt = indx[ir];
			q = arrin[indxt];
			indx[ir] = indx[0];
			if((--ir)<=0) {
				indx[0] = indxt;
				return ;
			}
		}
		i = L;
		j = L + L + 1;
		while (j<=ir) {
			if (j<ir)
				if(arrin[indx[j]] < arrin[indx[j+1]]) j++;
			if (q<arrin[indx[j]]) {
				indx[i] = indx[j];
				i = j;
				j = j + j + 1;
			}
			else j = ir + 1 ;
		}
		indx[i] = indxt;
	} // for(;;)
}
 


