#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double* pbootsamp = NULL;
    double* pp_y = NULL;
    double* po_y = NULL;
    double* pB = NULL;
    int bootRow = 0, pRow = 0;
    int bootCol = 0, pCol = 0;
    int i=0,j=0,index=0,m=0,n=0,t = 0, s = 0;
    double tem_sum = 0,ymean = 0, xmean = 0, tem_swap = 0;
    double ssx = 0;
    double tema = 0, temb = 0, suma = 0, sumb = 0;
    double x[9],y[9],tem_y[9],or_y[9];
    double es_tem, es_temp, es_sum;
    int es_up, es_down;
    if (nrhs != 3)
    {
        mexErrMsgTxt("Wrong inputs!!");
    }
    pbootsamp = mxGetPr(prhs[0]);                           //bootsamp
    pp_y = mxGetPr(prhs[1]);                                //log-picture data
    po_y = mxGetPr(prhs[2]);                                //original data
    bootRow = mxGetM(prhs[0]);                              //bootsamp row
    bootCol = mxGetN(prhs[0]);                              //bootsamp column
    pRow = mxGetM(prhs[1]);                                 //picture row              
    pCol = mxGetN(prhs[1]);                                 //picture column
    plhs[0] = mxCreateDoubleMatrix(pRow, pCol, mxREAL);     //create output variable
    pB = mxGetPr(plhs[0]);                                  //get the pointer of the output
    
    for(i = 1; i < (pCol-1); i++)
    {
        for(j = 1; j < (pRow-1); j++)
        {
            index = 0;
            suma = 0;
            sumb = 0;
            for(m = -1; m < 2; m++)
            {
                for(n = -1; n < 2; n++)
                {
                    tem_y[index] = pp_y[(i + m) * pRow + j + n];
                    or_y[index] = po_y[(i + m) * pRow + j + n];
                    index++;
                }
            }
            if(or_y[4] != 0)
            {
                t = 0;
                while(1)
                {
                    es_sum = 1;
                    es_temp = 1;
                    if (t != 0)
                    {
                        for(m = 1; m <= t; m++)
                        {
                            es_tem = 1;
                            for(n = 1; n <= m; n++)
                            {
                                es_tem = es_tem*or_y[4]/n;
                            }
                            es_sum = es_sum + es_tem;
                        }
                    }
                    else
                    {
                        es_sum = 1;
                    }
                    for(m = 1; m <= or_y[4]; m++)
                    {
                        es_temp = es_temp * 2.71828;
                    }
                    es_temp = es_sum / es_temp;
                    if(es_temp > 0.025)
                    {
                        es_down = t;
                        break;
                    }
                    t++;
                }
                t = 1;
                while(1)
                {
                    es_sum = 1;
                    es_temp = 1;
                    for(m = 1; m <= t; m++)
                    {
                        es_tem = 1;
                        for(n = 1; n <= m; n++)
                        {
                            es_tem = es_tem*or_y[4]/n;
                        }
                        es_sum = es_sum + es_tem;
                    }
                    for(m = 1; m <= or_y[4]; m++)
                    {
                        es_temp = es_temp * 2.71828;
                    }
                    es_temp = es_sum / es_temp;
                    if(es_temp > 0.975)
                    {
                        es_up = t;
                        break;
                    }
                    t++;
                }
                
                for(m = 0; m < 9; m++)
                {
                    if(or_y[m] < es_down)
                    {
                        tem_y[m] = tem_y[4];
                    }
                    if(or_y[m] > es_up)
                    {
                        tem_y[m] = tem_y[4];
                    }
                }
            }
            
            for(m = 0; m < 9; m++)
            {
                for(n = 0; n < (8 - m); n++)
                {
                    if(tem_y[n] > tem_y[n+1])
                    {
                        tem_sum = tem_y[n+1];
                        tem_y[n+1] = tem_y[n];
                        tem_y[n] = tem_sum;
              //          tem_sum = orx[n+1];
              //          orx[n+1] = orx[n];
              //          orx[n] = tem_sum;
                    }
                }
            }
            for(m = 0; m < bootCol; m++)
            {
                tem_sum = 0;
                for(n = 0; n < bootRow; n++)
                {
                    index = pbootsamp[m*bootRow+n] - 1;
                    y[n] = tem_y[index];
                    x[n] = index + 1;
                    tem_sum = tem_sum + x[n]*y[n];
                }
                /*for(t = 0; t < 9; t++)
                {
                    for(s = 0; s < (8 - t); s++)
                    {
                        if(y[s] > y[s+1])
                        {
                            tem_swap = y[s+1];
                            y[s+1] = y[s];
                            y[s] = tem_swap;
                        }
                    }
                }*/
               // for(n = 0; n < bootRow; n++)
               // {
               //     tem_sum = tem_sum + x[n]*y[n];
               // }
                xmean = (x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8])/9;
                ymean = (y[0]+y[1]+y[2]+y[3]+y[4]+y[5]+y[6]+y[7]+y[8])/9;
                ssx = (x[0]-xmean)*(x[0]-xmean)+(x[1]-xmean)*(x[1]-xmean)+(x[2]-xmean)*(x[2]-xmean)+(x[3]-xmean)*(x[3]-xmean)+(x[4]-xmean)*(x[4]-xmean)+(x[5]-xmean)*(x[5]-xmean)+(x[6]-xmean)*(x[6]-xmean)+(x[7]-xmean)*(x[7]-xmean)+(x[8]-xmean)*(x[8]-xmean);
                tema = (tem_sum - 9*xmean*ymean)/ssx;
                temb = ymean - tema*xmean;
                suma = suma + tema;
                sumb = sumb + temb;
            }
            pB[i * pRow + j] = 5 * suma / bootCol + sumb / bootCol;
            /*for(n = 0; n < 9; n++)
            {
                if(orx[n] == 5)
                {
                    tem_swap = (n+1) * suma / 200 + sumb / 200;
                    if(tem_swap < 3.5607)
                    {
                        pB[i * pRow + j] = (n+1) * suma / 200 + sumb / 200;
                    }
                    else
                    {
                        pB[i * pRow + j] = pp_y[i * pRow + j];
                    }
                }
            }*/
        }
    }
    //pnum[0] = pnum1[0] + pnum2[0];
    //mexPrintf("Done!\n");
}