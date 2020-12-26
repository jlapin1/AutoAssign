/*
 * This version reads in precalculated crosspeak matrices
 * as well as consolidates the peak and object structures
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <dirent.h>

///////////////////////////
///Global Variables////////
///////////////////////////
int 
	steps = 1E5,
	refine = 0,
	top =100000,//Number of top scores saved at the end
	outp = 1,// 1 for extended output, 0 for final results only
	otol1 = 2,
	otol2 = 0,
	soltol = 1E5;
double
	Wmin = 1E-8,
	Wmax = 0.5E1,
	no_incr = 1E3,
	alpha = 0.95,
	ocut1 = 0.1,
	ocut2 = 0.1,
    w1 = 1,
    w2 = 1,
	no_incr2,
	W;

double
	E_1_hold = 0,
	E_1_new = 0,
	E_2_hold = 0,
	E_2_new = 0,
	E_hold = 0,
	E_new = 0,
	E_low = 0,
	E1_low,
	E2_low,
	diff,
	diff2;

char
	* filepath = "Potgrids/",
	* outpath = "Output/",
	//filename1[20],
	* full,
	* full2;

int i,j,k,l,m,n,o,
    	amt,
	accept_l,
	accept_h,
	reject,
	first,
	second,
	tick = 0,
	tick2 = 0,
	rtick = 0,
	ftick = 0,
	solutick = 0,
	flor,
	ceel;

char
	filename1[20] = "input.csv";

typedef struct{
	FILE * fp;
	char
		*res;
	int 
        amt,
		len,
        xlen,
		**pks;
	double
		***E;

} inpobj;

typedef struct{
    int 
        **pkscur, **pkshold,
        **xpks;
    double
        **rEcur, **rEhold,
        *weight,
        *Ecur, *Ehold;
} pkobj;

// Anchor Parameters
char
	anchorq1,
	anchorq2,
	** aa;

int
	* res_cnt;

clock_t start_program, end_program;
///////////////////////////////////
///End global variables////////////
///////////////////////////////////


/////////////////////////////
/// Functions/subroutines ///
/////////////////////////////
void copyRxC(int **dest, int **src, int rows, int cols){
    int x,y;
    for(x=0;x<rows;x++)
        for(y=0;y<cols;y++)
            dest[x][y] = src[x][y];
}
void readinput(char *fpath, inpobj * a, pkobj * b){
	int x,y,z;
	a->fp = fopen(fpath,"r");
	if(a->fp==NULL)
		perror("File not found: ");
    
    // amt
    first = fscanf(a->fp, "%d\n", &a->amt);
    
    // len
	first = fscanf(a->fp,"%d\n",&a->len);//works
    
    // pks
	a->pks = (int**)malloc(2*sizeof(int*));
    a->pks[0] = (int*)malloc(a->len*sizeof(int));
    first = fscanf(a->fp,"%d",&a->pks[0][0]);
    for(y=1;y<a->len;y++)
        first = fscanf(a->fp,",%d",&a->pks[0][y]);
    first = fscanf(a->fp,"\n");
	a->pks[1] = (int*)calloc(a->len,sizeof(int));
	
    // res
    a->res = (char*)malloc(a->len*sizeof(char));
	first = fscanf(a->fp,"%c",&a->res[0]);
	for(x=1;x<a->len;x++)
		first = fscanf(a->fp,",%c",&a->res[x]);
	first = fscanf(a->fp,"\n");
	
    // E
    a->E = (double***)malloc(a->amt*sizeof(double**));
    for(x=0; x<a->amt; x++){
        a->E[x] = (double**)malloc(a->len*sizeof(double*));
        for(y=0; y < a->len; y++){
            a->E[x][y] = (double*)calloc(a->len, sizeof(double));
            first = fscanf(a->fp, "%lf", &a->E[x][y][0]);
            for(z=1; z<a->len; z++)
                first = fscanf(a->fp, ",%lf", &a->E[x][y][z]);
            first = fscanf(a->fp, "\n");
        }
	}
	fclose(a->fp);
    
    // Allocate for the pkobj
    a->xlen = 2*(a->len-1);
    b->pkscur = (int**)malloc(2*sizeof(int*));
    b->pkshold = (int**)malloc(2*sizeof(int*));
    b->xpks = (int**)malloc(2*sizeof(int*));
    for(x=0;x<2;x++){
        b->pkscur[x] = (int*)calloc(a->len, sizeof(int));
        b->pkshold[x] = (int*)calloc(a->len, sizeof(int));
        b->xpks[x] = (int*)calloc(a->xlen, sizeof(int));
    }
    b->rEcur = (double**)malloc(a->amt*sizeof(double*));
    b->rEhold = (double**)malloc(a->amt*sizeof(double*));
    for(x=0; x<a->amt; x++){
        b->rEcur[x] = (double*)calloc(a->xlen, sizeof(double));
        b->rEhold[x] = (double*)calloc(a->xlen, sizeof(double));
    }
    b->Ecur = (double*)calloc(a->amt+1, sizeof(double));
    b->Ehold = (double*)calloc(a->amt+1, sizeof(double));
}
void readconf(char *fpath){
	// This subroutine sets a bunch of global variables
	// atoll(): convert string to integer
	// atof(): convert string to float
	FILE * u = fopen(fpath, "r");
	char * nm = (char*)malloc(20*sizeof(char)),
		 * opt = (char*)malloc(20*sizeof(char));
	int fout;
	while(feof(u)==0){
		fout = fscanf(u,"%s", nm); // Read first string
		if(*nm!='#'){ // If not a comment line
			fout = fscanf(u,"%s\n",opt); // Read second string
			if (strcmp(nm, "amt")==0)
				amt = atoll(opt);
			else if(strcmp(nm,"steps")==0)
				steps = (int)atof(opt); // BUG: atoll("1E5") --> 1; must read in as a float and cast to integer
			else if(strcmp(nm,"top")==0)
				top = (int)atof(opt);
			else if(strcmp(nm,"outp")==0)
				outp = atoll(opt);
			else if(strcmp(nm,"otol1")==0)
				otol1 = (int)atof(opt);
			else if(strcmp(nm,"otol2")==0)
				otol2 = atoll(opt);
			else if(strcmp(nm,"soltol")==0)
				soltol = (int)atof(opt);
			else if(strcmp(nm,"Wmin")==0)
				Wmin = atof(opt);
			else if(strcmp(nm,"Wmax")==0)
				Wmax = atof(opt);
			else if(strcmp(nm,"no_incr")==0)
				no_incr = atof(opt);
			else if(strcmp(nm,"alpha")==0)
				alpha = atof(opt);
			else if(strcmp(nm,"ocut1")==0)
				ocut1 = atof(opt);
			else if(strcmp(nm,"ocut2")==0)
				ocut2 = atof(opt);
            else if(strcmp(nm,"w1")==0)
                w1 = atof(opt);
            else if(strcmp(nm,"w2")==0)
                w2 = atof(opt);
		}
		else // If a comment (#), skip everything that is not a newline and then read the newline
			fout = fscanf(u,"%*[^\n]\n");
	}
	fclose(u);
	free(nm);free(opt);
}
void xpeaks(int length, int ** pks1, int ** xpks1){
	int x,y,z;
	
	for(x=0;x<(length-1);x++){
		xpks1[0][x] = pks1[0][x];
		xpks1[1][x] = pks1[0][x+1];

		xpks1[0][x+length-1] = pks1[0][x+1];
		xpks1[1][x+length-1] = pks1[0][x];
	}
}
void xnumbers(inpobj a, int * inds, int * num, int * xnums){
	/*
	This function, based on the index swap, will find the number of crosspeak indices that will change their crosspeak (num),
	and their crosspeak indices (xnums). The output can be used as input to eng3.
	*/
	//7 Scenarios
	if(inds[0]==0){//First Left end
		if(inds[1]==1){//Second Adjacent
			*num = 4;
			xnums[0]=0;xnums[1]=1;xnums[2]=0+a.len-1;xnums[3]=1+a.len-1;
		}
		else if(inds[1]==(a.len-1)){//Second Right end
			*num=4;
			xnums[0]=0;xnums[1]=a.len-2;xnums[2]=0+a.len-1;xnums[3]=(2*a.len)-3;
		}
		else{//Second Middle
			*num=6;
			xnums[0]=0;xnums[1]=inds[1]-1;xnums[2]=inds[1];xnums[3]=0+a.len-1;xnums[4]=inds[1]+a.len-2;xnums[5]=inds[1]+a.len-1;
		}
	}
	else if(inds[0]==(a.len-2)){//Adjacent to right end
		*num=4;
		xnums[0]=a.len-3;xnums[1]=a.len-2;xnums[2]=(2*a.len)-4;xnums[3]=(2*a.len)-3;
	}
	else{//First Middle
		if(abs(inds[0]-inds[1])==1){//Second Adjacent
			*num=6;
			xnums[0]=inds[0]-1;xnums[1]=inds[0];xnums[2]=inds[1];xnums[3]=inds[0]+a.len-2;xnums[4]=inds[0]+a.len-1;xnums[5]=inds[1]+a.len-1;
		}
		else if(inds[1]==(a.len-1)){//Second Right end
			*num=6;
			xnums[0]=inds[0]-1;xnums[1]=inds[0];xnums[2]=inds[1]-1;xnums[3]=inds[0]+a.len-2;xnums[4]=inds[0]+a.len-1;xnums[5]=inds[1]+a.len-2;
		}
		else{//Both in the middle and not adjacent
			*num=8;
			xnums[0]=inds[0]-1;xnums[1]=inds[0];xnums[2]=inds[1]-1;xnums[3]=inds[1];xnums[4]=inds[0]+a.len-2;xnums[5]=inds[0]+a.len-1;xnums[6]=inds[1]+a.len-2;xnums[7]=inds[1]+a.len-1;
		}
	}
}
double* linspace(double lower, double upper, int points){
	// Must deallocate memory if this function is called many times
	int x,y,z;
	double 
		spacing,
		* output = (double*)malloc(points*sizeof(double));

	spacing = (upper-lower) / (double)(points-1);

	for(x=0;x<points;x++)
		output[x] = lower + (spacing*x);

	return output;

}
void eng2(inpobj a, pkobj * b){
	// Potential speed up if you calculate half the crosspeaks and multiply by two
	int x,y;
    double w[2] = {w1, w2};
    
	for(x=0;x<a.amt;x++){
        b->Ecur[x] = 0;
        for(y=0; y<a.xlen; y++){
            b->rEcur[x][y] = a.E[x][b->xpks[1][y]][b->xpks[0][y]];
            b->Ecur[x] += b->rEcur[x][y];
        }
        b->Ecur[a.amt] += w[x]*b->Ecur[x];
	}
}
void eng3(inpobj a, pkobj * b, int * xnums, int num){
	// Speeds up energy calculation over eng2 by updating energy based on swap, instead of re-calculating sum of all crosspeaks
	int x,y;
    double w[2] = {w1,w2};

	//Must copy old crosspeak energies into new crosspeak energies (not too costly)
	for(x=0;x<a.amt;x++){
        b->Ecur[x] = b->Ehold[x];
        for(y=0;y<a.xlen;y++)
            b->rEcur[x][y] = b->rEhold[x][y];
	}
	b->Ecur[a.amt] = 0;
    for(x=0;x<a.amt;x++){
        for(y=0;y<num;y++){
	        b->Ecur[x] -= b->rEhold[x][xnums[y]];
            b->rEcur[x][xnums[y]] = a.E[x][b->xpks[1][xnums[y]]][b->xpks[0][xnums[y]]];
            b->Ecur[x] += b->rEcur[x][xnums[y]];
        }
        b->Ecur[a.amt] += w[x]*b->Ecur[x];
	}
}
void swap_ind(int length, int * res_count, int ** pks1_in, int ** pks1_out, int * inds){
	int x, y, z;
	int
		indo = 10000,
		indn = 10000,
		indnn,
		tick = 0,
		tick2 = 0,
		used,
        hold1,
        hold2,
        hold3,
        hold4;
    
    copyRxC(pks1_out, pks1_in, 2, length);

	indo = rand()%length;//Seed indo
	while( (res_count[(int)pks1_out[1][indo]] == 1) || ((int)pks1_out[1][indo]==21) ) // If an index's residue only occurs once or it's a 'Z', that index will always have the same chemical shift
		indo = rand()%length;
	indn = rand()%length;//Seed indn
	while( (indo==indn) || (pks1_out[1][indo]!=pks1_out[1][indn]) || (pks1_out[1][indn]==21) )
		indn = rand()%length;		
	//Should have good indices for the swap at this point
	
	hold1 = pks1_out[0][indn];
	pks1_out[0][indn] = pks1_out[0][indo];
	pks1_out[0][indo] = hold1;
	
	if(indo<indn){
		inds[0] = indo;
		inds[1] = indn;
	}
	else{
		inds[0] = indn;
		inds[1] = indo;
	}
}
int* res_type(int length, char * res_inp, int ** pks1){
	int x, y, z;

	int * res_count = (int*)calloc(22,sizeof(int));
	
	for(x=0;x<length;x++){
		switch(res_inp[x]){
			case 'A'://Alanine
				pks1[1][x] = 1;
				res_count[1]++;
				break;
			case 'C'://Cysteine
				pks1[1][x] = 2;
				res_count[2]++;
				break;
			case 'D'://Aspartic acid
				pks1[1][x] = 3;
				res_count[3]++;
				break;
			case 'E'://Glutamine
				pks1[1][x] = 4;
				res_count[4]++;
				break;
			case 'F'://Phenylalanine
				pks1[1][x] = 5;
				res_count[5]++;
				break;
			case 'G'://Glycine
				pks1[1][x] = 6;
				res_count[6]++;
				break;
			case 'H'://Histidine
				pks1[1][x] = 7;
				res_count[7]++;
				break;
			case 'I'://Isoleucine
				pks1[1][x] = 8;
				res_count[8]++;
				break;
			case 'K'://Lysine
				pks1[1][x] = 9;
				res_count[9]++;
				break;
			case 'L'://Leucine
				pks1[1][x] = 10;
				res_count[10]++;
				break;
			case 'M'://Methionine
				pks1[1][x] = 11;
				res_count[11]++;
				break;
			case 'N'://Arginine
				pks1[1][x] = 12;
				res_count[12]++;
				break;
			case 'P'://Proline
				pks1[1][x] = 13;
				res_count[13]++;
				break;
			case 'Q'://Glutamine
				pks1[1][x] = 14;
				res_count[14]++;
				break;
			case 'R'://Arginine
				pks1[1][x] = 15;
				res_count[15]++;
				break;
			case 'S'://Serine
				pks1[1][x] = 16;
				res_count[16]++;
				break;
			case 'T'://Threonine
				pks1[1][x] = 17;
				res_count[17]++;
				break;
			case 'V'://Valine
				pks1[1][x] = 18;
				res_count[18]++;
				break;
			case 'W'://Tryptophan
				pks1[1][x] = 19;
				res_count[19]++;
				break;
			case 'Y'://Tyrosine
				pks1[1][x] = 20;
				res_count[20]++;
				break;
			case 'Z'://Fixed residue/shift
				pks1[1][x] = 21;
				res_count[21]++;
				break;
			default://X or other
				pks1[1][x] = 0;
				res_count[0]++;
		}
	}
	return res_count;
}
int orfilt(inpobj a, double * engs, double cut, int tol){
	int i, tal;
	tal=0;
	for(i=0;i<2*(a.len-1);i++)
		if(engs[i]>cut){
			tal++;
			if(tal>tol){
				return 0;
				break;
			}
		}
	if(tal<=tol)
		return 1;
	else
		return 0;
}
int xpkfilt(inpobj a, double **xpks, double ymin, double ymax, double xmin, double xmax, int num, int typ){
	int i, tally;
	tally = 0;
	for(i=0;i<2*(a.len-1);i++)
	{
		if( (ymin<=xpks[0][i]) && (ymax>=xpks[0][i]) && (xmin<=xpks[1][i]) && (xmax>=xpks[1][i]) )
			tally++;
		if(typ==-1) // At least
		{
			if(tally==num)
				return 1;//Quick exit
		}
		else if(typ==1) // At most
		{
			if(tally>num)
				return 0;//Quick exit
		}
	}
	if(typ==-1)
		return 0;
	if(typ==0)
	{
		if(tally==num)
			return 1;
		else
			return 0;
	}
	if(typ==1)
		return 1;
}
/////////////////////////////
///End of functions//////////
/////////////////////////////



int main(){ 
	// I will be using native C random number methods
	srand(time(NULL));
	
	// Look for configuration script "config.conf"
    struct dirent * de;
    DIR * dr = opendir("./");
    if(dr==NULL){
        printf("\nCould not open directory.\n");
        return 0;
    }
    int filexist = 0;
    while((de = readdir(dr)) != NULL){
        filexist = strcmp(de->d_name, "config.conf");
        if(filexist==0)
            break;
    }
    closedir(dr);
    
    // Read into configuration script to assign global variables
    if(filexist==1){
        printf("\nNo configuration file \"config.conf\" found in current directory\n\n");
        exit(0);
    }
    else
        readconf("config.conf");
	
	start_program = clock();
	///////////////////////////////////////////////////////////////////////////////////
	///Read in experimental pseudo-potental////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	
	FILE * f, * g, * h;
	
	inpobj
		ob1,
		ob2;
    
    pkobj
        pk1,
        pk2;

	first = strlen(filepath);
	second = strlen(filename1);
	full = (char*)malloc((first+second)*sizeof(char));
	strcpy(full,filepath);
	strcat(full,filename1);
	readinput(full, &ob1, &pk1);
	///////////////////////////////////////////////////////////////////////////////////
	///BREAK///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////////////
	/// Monte Carlo/Simulated Annealing ///////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	double start = log(1-(Wmin/Wmax))/log(alpha); // Exponential
	//double start = sqrt( ((1.0/(1-(Wmin/Wmax)))-1) / alpha ); // Quadratic
	//double start = Wmin; // Linear
	double * incr = linspace(start, 100, no_incr);
	if(refine!=0)
		no_incr2 = no_incr + ((refine) * (int)(0.9*no_incr));
	else
		no_incr2 = no_incr;

	int 
		* inds = (int*)malloc(2*sizeof(int)),
		* xnums = (int*)malloc(8*sizeof(int)),
		num;

	int
		** rpks_save = (int**)malloc(top*sizeof(int*));
    double
        ** E_save = (double**)malloc(top*sizeof(double*)),
		** E_save_all = (double**)malloc(4*sizeof(double*)),
		
		cur,
		acc_crit,
		dE;
	
	// Allocate for saving variables
	for(i=0;i<top;i++){
		rpks_save[i] = (int*)malloc(ob1.len*sizeof(int));
        E_save[i] = (double*)calloc(3, sizeof(double));
    }
    
	for(i=0;i<4;i++)
		E_save_all[i] = (double*)malloc((no_incr2+1)*sizeof(double)); //FIX
	
	// Start the fun
	// Copy peaks into other variables
    copyRxC(pk1.pkscur, ob1.pks, 2, ob1.len);
    
    // Assign residues to every peak
    // After this function is run, residue assignments to indices should be permanent throughout the program
    res_cnt = res_type(ob1.len, ob1.res, pk1.pkscur);
    
    // Shuffle the deck
    for(i=0;i<10000;i++)
        swap_ind(ob1.len, res_cnt, pk1.pkscur, pk1.pkscur, inds);
    
    // Generate crosspeaks
    xpeaks(ob1.len, pk1.pkscur, pk1.xpks);

    // Calculate pseudopotential of crosspeaks on E surface
    eng2(ob1, &pk1);
    
    copyRxC(pk1.pkshold, pk1.pkscur, 2, ob1.len);
    for(o=0;o<ob1.amt;o++){
        pk1.Ehold[o]=pk1.Ecur[o];
        for(i=0;i<ob1.xlen;i++)
            pk1.rEhold[o][i] = pk1.rEcur[o][i];
    }
    pk1.Ehold[ob1.amt] = pk1.Ecur[ob1.amt];
    E_hold = pk1.Ehold[ob1.amt];
    
    W = Wmin;
    E_save_all[0][0] = W;
    E_save_all[1][0] = pk1.Ecur[0];
    E_save_all[2][0] = pk1.Ecur[1];
    E_save_all[3][0] = pk1.Ecur[2];
    
    o=0;
    while(o<no_incr){
        accept_h = 0;
        accept_l = 0;
        reject = 0;
        
        cur = 100*((double)o/(double)no_incr);
        
        if(outp==1){
            if(refine>0)
                printf("Refine: %d\n",tick);
            printf("W = %.3lf (%.0lf%% complete)\n",W,cur);
        }

        solutick=0;
        //Take random step
        for(i=0;i<steps;i++){
            rtick = 0;
            ftick = 0;
            
            // New peak order
            //	- rpks will be the new try
            swap_ind(ob1.len, res_cnt, pk1.pkshold, pk1.pkscur, inds);
            
            //Setting num and xnums
            xnumbers(ob1, inds, &num, xnums);
            
            // Generate crosspeaks
            xpeaks(ob1.len, pk1.pkscur, pk1.xpks);
            
            // Compute pseudopotential
            eng3(ob1, &pk1, xnums, num);
            E_new = pk1.Ecur[ob1.amt];
            
            // E_low is for stdout only
            if(E_new<E_low){
                E_low=E_new;
                E1_low=pk1.Ecur[0];
                E2_low=pk1.Ecur[1];
            }
            
            if( E_new < E_hold ){ // Energy went down
                E_hold = E_new;
                copyRxC(pk1.pkshold, pk1.pkscur, 2, ob1.len);
                for(j=0;j<ob1.amt;j++){
                    pk1.Ehold[j] = pk1.Ecur[j];
                    for(k=0;k<ob1.xlen;k++)
                        pk1.rEhold[j][k]=pk1.rEcur[j][k];
                }
                pk1.Ehold[ob1.amt] = pk1.Ecur[ob1.amt];
                accept_l++;
            }
            else if(E_new == E_hold){
                acc_crit = (double)rand() / (double)RAND_MAX;
                if(acc_crit>0.5){//Accept equal energy momve half the time --> stdout tallies may not sum to steps --> expect < 1% to be equal
                    E_hold = E_new;
                    copyRxC(pk1.pkshold, pk1.pkscur, 2, ob1.len);
                    for(j=0;j<ob1.amt;j++){
                        pk1.Ehold[j] = pk1.Ecur[j];
                        for(k=0;k<ob1.xlen;k++)
                            pk1.rEhold[j][k]=pk1.rEcur[j][k];
                    }
                    pk1.Ehold[ob1.amt] = pk1.Ecur[ob1.amt];
                }
                else{
                    reject++;
                    rtick = 1;
                }
            }
            else{ // Energy went up  --> apply Metropolis criterion
                acc_crit = (double)rand() / (double)RAND_MAX;
                dE = E_new - E_hold;
                if( exp(-1*dE*W) > acc_crit ){
                    E_hold = E_new;
                    copyRxC(pk1.pkshold, pk1.pkscur, 2, ob1.len);
                    for(j=0;j<ob1.amt;j++){
                        pk1.Ehold[j] = pk1.Ecur[j];
                        for(k=0;k<ob1.xlen;k++)
                            pk1.rEhold[j][k]=pk1.rEcur[j][k];
                    }
                    pk1.Ehold[ob1.amt] = pk1.Ecur[ob1.amt];
                    accept_h++;
                }
                else{
                    reject++;
                    rtick = 1;
                }
            }//Random step is over

            
            // Sort Energy saving variables
            //	- Issue: same equal numbers are not registered by the == logical operator
            //		- The result of that issue caused reverse orders to be registered adjacent
            //		- Used a workaround by calculating differences between E_hold and E_save[m]				
            
            // Apply filter
            if( (orfilt(ob1, pk1.rEcur[0], ocut1, otol1)==1) && (orfilt(ob1, pk1.rEcur[1], ocut2, otol2)==1) )
                ftick=1;
            
            diff = 1;
            m=-1;
            flor = 0;
            ceel = top-1;
            j=round((ceel)/2);
            if( (E_new<E_save[ceel][ob1.amt]) && (ftick==1) ){
                if(j==0)
                    m=j;
                while( (j!=0) && (j!=(top-1)) ){//(E_hold<=E_save[j]) && (j>-1) && (rtick==0) )
                    if(E_new<E_save[j][ob1.amt]){
                        diff = pow(pow(E_new-E_save[j][ob1.amt],2.0),0.5); //If an "equal" number sneaks in the backdoor
                        diff2 = pow(pow(E_new-E_save[j-1][ob1.amt],2.0),0.5);
                        if( (E_new==E_save[j-1][ob1.amt]) || (diff<1E-6) || (diff2<1E-6)){
                            m=-1;
                            break;
                        }
                        if( (E_new>E_save[j-1][ob1.amt]) || (j==0) ){
                            m=j;
                            break;
                        }
                        ceel = j;
                        j = round((j+flor)/2);
                        if(j==0)
                            m=j;
                    }
                    else if(E_new>E_save[j][ob1.amt]){
                        diff = pow(pow(E_new-E_save[j][ob1.amt],2.0),0.5); //If an "equal" number sneaks in the backdoor
                        diff2 = pow(pow(E_new-E_save[j+1][ob1.amt],2.0),0.5);
                        if( (E_new==E_save[j][ob1.amt]) || (diff<1E-6) || (diff2<1E-6) ){
                            m=-1;
                            break;
                        }
                        if( E_new<E_save[j+1][ob1.amt] ){
                            m=j+1;
                            break;
                        }
                        flor = j;
                        j = round((j+ceel)/2);
                    }
                    else{
                        m=-1;
                        break;
                    }
                }
            }
            if(m!=-1){
                solutick++;
                for(n=(top-1);n>m;n--){ // Move everything at that index downstream 1 spot
                    for(j=0;j<ob1.amt+1;j++)
                        E_save[n][j] = E_save[n-1][j];
                    for(j=0;j<ob1.len;j++){
                        rpks_save[n][j] = rpks_save[n-1][j];
                    }
                }
                // Now insert the new value at index m
                for(j=0;j<ob1.amt+1;j++)
                    E_save[m][j] = pk1.Ecur[j];
                for(n=0;n<ob1.len;n++){rpks_save[m][n]=pk1.pkscur[0][n];}
                
            }// End saving energies/peaks
        }// W run is over

        if(outp==1){
            printf("\n\n");
            printf("	Steps accepted (E up)  : %d\n",accept_h);
            printf("	Steps accepted (E down): %d\n",accept_l);
            printf("	Steps rejected         : %d\n",reject);
            printf("	Current E: %lf\n",E_hold);
            if(amt==1)
                printf("	Low E	 : %lf\n",E_low);
            else if(amt>1)
                printf("	Low E    : %lf (hom = %lf, het = %lf)\n",E_low,E1_low,E2_low);
            j=0;
            while( (j<top) && (E_save[j][ob1.amt]!=0.0) )
                j++;
            printf("	Solutions: %d\n",j);
            printf("\n\n");
        }

        // Store information for energy profile
        E_save_all[0][o+1] = W;
        E_save_all[1][o+1] = pk1.Ehold[0];
        E_save_all[2][o+1] = pk1.Ehold[1];
        E_save_all[3][o+1] = pk1.Ehold[2];

        // Post cycle refinement by increasing temperature (W)
        if((o==no_incr-1) && (tick!=refine) ){
            o = (int)(no_incr/10);
            tick++;
        }

        W = Wmax*(1 - pow(alpha,incr[o]));
        //W = Wmax*(1 - (1/(1+(alpha*incr[o]*incr[o])))); //Quadratic
        //W = incr[o] //Linear
        
        if(solutick<soltol) // If less than soltol solutions were tallied, continue on to next temperature
            o++;
    }
	
	////////////////////////
	// Print results to file
	////////////////////////
	char outfile[15] = "output.csv";
	char outfull[100];
	strcpy(outfull,outpath);
	strcat(outfull,outfile);
	k=0;
	while((k<top) && (E_save[k][ob1.amt]!=0) )
		k++;
	g = fopen(outfull,"w");
	fprintf(g,"%d\n",k);
	for(m=0;m<k;m++){	
        fprintf(g, "%d",rpks_save[m][0]);
        for(j=1;j<ob1.len;j++)
            fprintf(g,",%d",rpks_save[m][j]);
        fprintf(g,"\n");
	}
    fclose(g);
	
	strcpy(outfull,outpath);
	strcat(outfull,"E_profile.csv");
	g = fopen(outfull,"w");
	fprintf(g,"%d\n",(int)no_incr);
	for(i=0;i<=no_incr;i++)
		fprintf(g,"%lf,%lf,%lf,%lf\n",E_save_all[0][i],E_save_all[1][i],E_save_all[2][i],E_save_all[3][i]);
	fclose(g);

	printf("\nTop Scores\n");
	for(i=0;i<k;i++)
		printf("%4d: %lf\n",i+1,E_save[i][ob1.amt]);
	
	end_program = clock();
	double elap = ((double)end_program-(double)start_program)/CLOCKS_PER_SEC;
	printf("\nTotal runtime: %.2lf seconds\n\n",elap);
    
}
