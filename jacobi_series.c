/* Jacobi迭代的串行实现 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
double A[1010][1010];
double B[1010];
double X[1010];
double X2[1010];
int main(){
    int n;
    printf("请输入矩阵阶数：\n");
    while(~scanf("%d",&n))
    {
        for(int i=1; i<=n; i++)    //输入A矩阵；
        {
            printf("请输入矩阵第%d行：\n",i);
            for(int j=1; j<=n; j++)
            {
                scanf("%lf",&A[i][j]);
            }
        }
 
        printf("请输入向量b：\n");
        for(int i=1; i<=n; i++) //输入B向量；
        {
            scanf("%lf",&B[i]);
        }
        memset(X,0,sizeof(X));
        memset(X2,0,sizeof(X));
 
        for(int i=0; i<13; i++) //雅可比迭代；
        {
            for(int j=1; j<=n; j++)
            {
                for(int k=1; k<=n; k++)
                {
                    if(j!=k)
                    {
                        X[j]+=(-1)*(A[j][k]/A[j][j])*X2[k];
                    }
                }
                X[j]+=B[j]/A[j][j];
 
            }
            printf("第%d次迭代结果\n",i+1);
            for(int j=1; j<=n; j++)
            {
                X2[j]=X[j];
                X[j]=0;
                printf("%.4f ",X2[j]);
            }
            printf("\n");
        }
        printf("请输入下一个矩阵阶数：\n");
    }
 
    return 0;
}
