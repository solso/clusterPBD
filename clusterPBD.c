#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <math.h>
#include <sys/timeb.h>
#include <pthread.h>

int** graph;
double** graph_weight;

int* conn;
int* self;
int graphN;
int graphE;
int P;
int D;
int mSave;
int mNumThreads;

int* partitions_to_save_v;
int partitions_to_save;

typedef struct {
    int i;
    int j;
    double eij;
    void* next;
    void* prev;
} edgeP;


int* sortDesc(int* keys, int* values, int size) {

    int* newkeys = (int *)malloc(sizeof(int)*size);
    int* taken = (int *)malloc(sizeof(int)*size);
    int i, j;
    int pos;
    int max;

    for(i=0;i<size;i++) taken[i] = 0;

    for(i=0;i<size;i++) {

        max = -999999999;
        pos = -1;
        for(j=0;j<size;j++) {
            if ((!taken[j]) && (values[j]>=max)) {
                max = values[j];
                pos = j;
            }
        }

        taken[pos] = 1;
        newkeys[i] = keys[pos];
    }

    free(taken);

    return newkeys;

}


void displayG() {
    int i,j;
    for(j=0;j<graphN;j++) {
        for(i=0;i<conn[j];i++) {
            printf("(%d,%d)\n",j,graph[j][i]);
        }
    }
}

void display(double** MAT, int n, int p) {
    int i,j;
    for(j=0;j<p;j++) {
        for(i=0;i<n;i++) {
            printf("%.3f ",MAT[i][j]);
        }
        printf("\n");
    }
}

void displayV(double* V, int n) {
    int j;
    for(j=0;j<n;j++) {
        printf("%.3f ",(double)V[j]);
    }
    printf("\n");
}


double modularity(int* part, int n) {
    int i;
    int j;
    int max;
    max=-1;

    for(i=0;i<n;i++) if (part[i]>max) max=part[i];
    max=max+1;

    int numlinks=0;

    int* ar = (int *)malloc(sizeof(int)*max);
    int* er = (int *)malloc(sizeof(int)*max);
    for(i=0;i<max;i++) {
        ar[i]=0;
        er[i]=0;
    }
    double q = 0.0;
    int id1, id2;

    for(i=0;i<n;i++) {
        id1=i;
        for(j=0;j<conn[i];j++) {
            id2=graph[i][j];
            if (part[id1]==part[id2]) er[part[id1]]++;
            ar[part[id1]]++;
            numlinks++;
            if (id1==id2) {
                ar[part[id1]]++;
                er[part[id1]]++;
                numlinks++;
            }
        }
    }

    for(i=0;i<max;i++) q = q + (((double)er[i]/(double)numlinks)-pow(((double)ar[i]/(double)numlinks),2.0));

    //printf("%d numlinks\n",numlinks);

    free(ar);
    free(er);

    return q;
}




int* distances(int id) {
    int i,j;
    int* dist = (int *)malloc(sizeof(int)*graphN);
    int* trav = (int *)malloc(sizeof(int)*graphN*2);
    for(i=0;i<graphN;i++) {
        dist[i]=0;
        trav[i]=-1;
    }

    int pos=0;
    int totEl=1;
    trav[pos]=id;

    while(pos<totEl) {
        i = trav[pos];
        pos++;
        for(j=0;j<conn[i];j++) {
            if ((dist[graph[i][j]]==0) && (graph[i][j]!=id)) {
                trav[totEl]=graph[i][j];
                totEl++;
                dist[graph[i][j]]=dist[i]+1;
            }
        }
    }

    for(i=0;i<graphN;i++) {
        if (dist[i]==0) {
            if (i!=id) {
                dist[i]=graphN*10;
            }
        }
    }

    free(trav);

    return dist;

}

int* neighbours(int ag1, int distance, int* size) {
    int i,j;
    int* dist = (int *)malloc(sizeof(int)*graphN);
    int* trav = (int *)malloc(sizeof(int)*graphN*2);

    for(i=0;i<graphN;i++) {
        trav[i]=-1;
        dist[i]=0;
    }

    int pos=0;
    int totEl=1;
    trav[pos]=ag1;

    while(pos<totEl) {
        i = trav[pos];
        if (dist[i]<distance) {
            //spos++;
            for(j=0;j<conn[i];j++) {
                if ((dist[graph[i][j]]==0) && (graph[i][j]!=ag1)) {
                    trav[totEl]=graph[i][j];
                    totEl++;
                    dist[graph[i][j]]=dist[i]+1;
                }
            }
        }
        pos++;
    }

    int cont=0;
    for(i=0;i<graphN;i++) {
        if (dist[i]!=0) cont++;
    }

    int* neig = (int *)malloc(sizeof(int)*cont);
    (*size)=cont;
    j=0;
    for(i=0;i<graphN;i++) {
        if (dist[i]!=0) {
            neig[j]=i;
            j++;
        }
        cont++;
    }

    free(dist);
    free(trav);

    return neig;

}

int* connexComp() {

    int i,j;
    int* dist = NULL;
    int* group = (int *)malloc(sizeof(int)*graphN);

    for(i=0;i<graphN;i++) group[i]=-1;

    int currGroup = 1;
    for(i=0;i<graphN;i++) {
        if (group[i]==-1) {
            dist = distances(i);
            for(j=0;j<graphN;j++) {
                if (dist[j]<graphN) {
                    group[j]=currGroup;
                }
            }
            free(dist);
            currGroup++;
        }
    }

    printf("connexComp: %d\n",currGroup-1);

    return group;

}

int* histo(int* v, int n, int* num) {
    int i;
    int max=-1;

    for(i=0;i<n;i++) {
        if (max<v[i]) max=v[i];
    }
    max=max+1;
    int* hist = (int *)malloc(sizeof(int)*max);

    for(i=0;i<max;i++) hist[i]=0;

    for(i=0;i<n;i++) hist[v[i]]++;

    (*num)=max;
    return hist;

}


int* matrixSeedMaxConn(int n, int* p) {


    int max;
    int* hist = histo(conn,n,&max);
    int i, j;
    int connthres=-1;
    int thres = (n / 5)+1;

    int cont;
    int acum=0;

    for(i=max-1;i>=0;i--) {
        //printf("acum %d, hist[%d]=%d\n",acum,i,hist[i]);
        if (acum<=thres) acum+=hist[i];
        else {
            connthres = i;
            break;
        }
    }

    if (connthres<=1) connthres=2;

    cont=0;
    for(i=0;i<n;i++) if (conn[i] >= connthres) cont++;

    int* res = (int *)malloc(sizeof(int)*cont);

    j=0;
    for(i=0;i<n;i++) {
        if (conn[i] >= connthres) {
            res[j]=i;
            j++;
        }
    }


    (*p)=cont;
    free(hist);

    return res;

}

int* load_initial_partition(char* filename, int* P) {

    FILE* file = fopen(filename,"r");
    char line[255];
    int id;
    char* s;
    int n = 0;
    int i = 0;
    int* res;
    int maxid = -1;

    while (!feof(file)) {
        s = fgets(line,255,file);
        n+=1;
    }
    n--;

    rewind(file);

    if (n!=graphN) {
        printf("Fatal error, n is not graphN, %d %d\n",n,graphN);
        exit(-1);
    }

    res = (int *)malloc(sizeof(int)*n);

    while (!feof(file)) {
        s = fgets(line,255,file);
        sscanf(line,"%d\n",&id);
        res[i]=id;
        if (id>maxid) maxid=id;
        i+=1;
    }

    fclose(file);

    (*P)=maxid+1;
    return res;

}

void load(char* filename) {

    FILE* file = fopen(filename,"r");
    char line[255];
    char temp[255];
    int id1;
    int id2;
    double w;
    int i, j;
    fpos_t pos;
    int n;
    char* s;

    graphE=0;

    s = fgets(line,255,file);
    sscanf(line,"%s %d\n",temp,&n);

    for(i=0;i<n;i++) {
        s = fgets(line,255,file);
    }

    s = fgets(line,255,file);


    conn = (int *)malloc(sizeof(int)*n);
    self = (int *)malloc(sizeof(int)*n);
    int* punt = (int *)malloc(sizeof(int)*n);

    for(i=0;i<n;i++) {
        conn[i]=0;
        punt[i]=0;
        self[i]=0;
    }

    fgetpos(file,&pos);
    int tt=0;
    while (!feof(file)) {
        s = fgets(line,255,file);
        if (strlen(line)>0) {
            sscanf(line,"%d %d %lf\n",&id1,&id2,&w);
            id1--;
            id2--;
            conn[id1]++;
            conn[id2]++;
            tt=tt+2;
        }
        line[0]='\0';
    }

    fsetpos(file,&pos);

    graph = (int **)malloc(sizeof(int *)*n);
    graph_weight = (double **)malloc(sizeof(double *)*n);
    for(i=0;i<n;i++) {
        graph[i] = (int *)malloc(sizeof(int)*(conn[i]));
        graph_weight[i] = (double *)malloc(sizeof(double)*(conn[i]));
        for(j=0;j<conn[i];j++) graph[i][j]=-1;
        for(j=0;j<conn[i];j++) graph_weight[i][j]=0.0;
    }

    int contli=0;
    int found;
    while (!feof(file)) {
        s = fgets(line,255,file);

        if (strlen(line)>0) {
            sscanf(line,"%d %d %lf\n",&id1,&id2,&w);
            id1--;
            id2--;


            found=0;
            for(i=0;i<punt[id1];i++) {
                if (graph[id1][i]==id2) found=1;
            }
            for(i=0;i<punt[id2];i++) {
                if (graph[id2][i]==id1) found=1;
            }

            if (found==0) {
                if (id1!=id2) {
                    graph[id1][punt[id1]]=id2;
                    graph_weight[id1][punt[id1]]=w;
                    punt[id1]++;

                    graph[id2][punt[id2]]=id1;
                    graph_weight[id2][punt[id2]]=w;
                    punt[id2]++;

                    //tt=tt+2;
                    graphE=graphE+1;
                }
                else {
                    //graph[id1][punt[id1]]=id2;
                    //punt[id1]++;
                    //graphE=graphE+1;
                    //self[id1]=1;
                }
            }
        }
        line[0]='\0';
        contli++;
    }

    graphN=n;

    fclose(file);

    for(i=0;i<n;i++) {
        conn[i]=punt[i];

    }

    int j2,k;
    for(i=0;i<n;i++) {
        for(j=0;j<conn[i];j++) for(j2=j+1;j2<conn[i];j2++) if (graph[i][j]==graph[i][j2]) {
            printf("node %d\n",i);
            for(k=0;k<conn[i];k++) printf(" %d",graph[i][k]);
            printf("\n");
            printf("eoeoeooe   \n"); exit(-1);
        }
    }

    free(punt);


}


void loadShort(char *filename) {

    FILE* file = fopen(filename,"r");

    char line[10000];
    char line2[10000];
    char temp[10000];
    int id1=-1;
    int id2=-1;
    int cont = 0;
    int i, j;
    int n;
    char* s;

    s = fgets(line,10000,file);
    sscanf(line,"%s %d\n",temp,&n);
    s = fgets(line,10000,file);

    graphN=n;
    graphE=0;

    conn = (int *)malloc(sizeof(int)*n);
    self = (int *)malloc(sizeof(int)*n);
    int* punt = (int *)malloc(sizeof(int)*n);

    for(i=0;i<n;i++) {
        conn[i]=0;
        punt[i]=0;
        self[i]=1;
    }

    while (!feof(file)) {
        s = fgets(line,10000,file);
        if (strlen(line)>0) {
            int first=1;
            int i = 0;
            int final = 0;
            line2[0] = '\0';

            while((!final) && (i<(int)strlen(line))) {

                if (((char)line[i] >= '0') && ((char)line[i] <= '9'))  {
                    sprintf(line2,"%s%c",line2,line[i]);
                }
                else if ((line[i] == ' ') || (line[i]== '\n')) {
                    if (strlen(line2)>0) {
                        //printf(":%s %d\n",line2,first);
                        if (first) {
                            id1 = atoi(line2);
                            if (id1==graphN) id1=0;
                            first=0;
                        }
                        else  {
                            id2 = atoi(line2);
                            if (id2==graphN) id2=0;

                            //printf("(%d,%d) ",id1,id2);
                            if (id1==id2) {
                                //printf("self (%d,%d)\n",id1,id2);
                            }
                            conn[id1]++;
                            conn[id2]++;
                        }
                        line2[0] = '\0';
                    }
                }
                else final = 1;
                i++;
            }

            line2[0] = '\0';
            //printf("\n");

            cont++;
        }

    }

    rewind(file);
    cont = 1;

    s = fgets(line,10000,file);
    s = fgets(line,10000,file);

    graph = (int **)malloc(sizeof(int *)*n);
    graph_weight = (double **)malloc(sizeof(double *)*n);

    for(i=0;i<n;i++) {
        graph[i] = (int *)malloc(sizeof(int)*(conn[i]));
        graph_weight[i] = (double *)malloc(sizeof(double)*(conn[i]));
        for(j=0;j<conn[i];j++) graph[i][j]=-1;
        for(j=0;j<conn[i];j++) graph_weight[i][j]=0.0;
    }

    int found, ik1;
    graphE=0;

    while (!feof(file)) {

        s = fgets(line,10000,file);

        if (strlen(line)>0) {

            int i = 0;
            int final = 0;
            int first = 1;

            line2[0] = '\0';

            while((!final) && (i<(int)strlen(line))) {

                if (((char)line[i] >= '0') && ((char)line[i] <= '9'))  {
                    sprintf(line2,"%s%c",line2,line[i]);
                }
                else if ((line[i] == ' ') || (line[i]== '\n')) {
                    if (strlen(line2)>0) {

                        if (first) {
                            id1 = atoi(line2);
                            if (id1==graphN) id1=0;
                            first=0;
                        }
                        else {
                            id2 = atoi(line2);
                            if (id2==graphN) id2=0;

                            if (id1==id2) {
                                //printf("self (%d,%d)\n",id1,id2);
                            }


                            found=0;
                            for(ik1=0;ik1<punt[id1];ik1++) {
                                if (graph[id1][ik1]==id2) found=1;
                            }
                            for(ik1=0;ik1<punt[id2];ik1++) {
                                if (graph[id2][ik1]==id1) found=1;
                            }

                            if (found==0) {
                                if (id1!=id2) {
                                    graph[id1][punt[id1]]=id2;
                                    punt[id1]++;
                                    graph[id2][punt[id2]]=id1;
                                    punt[id2]++;
                                    //tt=tt+2;
                                    graphE=graphE+1;
                                }
                                else {
                                    //graph[id1][punt[id1]]=id2;
                                    //punt[id1]++;
                                    //graphE=graphE+1;
                                    //self[id1]=1;
                                }
                            }


                        }

                        line2[0] = '\0';

                    }

                }
                else final = 1;

                i++;

            }



            line2[0] = '\0';




        cont++;



        }



   }

    for(i=0;i<n;i++) {
        conn[i]=punt[i];
    }


    int j2,k;
    for(i=0;i<n;i++) {
        for(j=0;j<conn[i];j++) for(j2=j+1;j2<conn[i];j2++) if ((graph[i][j]==graph[i][j2]) || (i==graph[i][j])) {
            printf("node %d\n",i);
            for(k=0;k<conn[i];k++) printf(" %d",graph[i][k]);
            printf("\n");
            printf("eoeoeooe   \n"); exit(-1);
        }
    }


    free(punt);

    fclose(file);



}



double norm(double* v, int n) {
    int i;
    double nor=0.0;
    for(i=0;i<n;i++) {
            nor+=v[i]*v[i];
    }
    nor=pow(nor,0.5);
    return nor;

}

double distance(double* v1, double* v2, int n) {
    int i;
    double nor=0.0;
    for(i=0;i<n;i++) {
        nor+=v1[i]*v2[i];
    }
    nor=1.0-(nor/(norm(v1,n)*norm(v2,n)));

    return nor;

}

double distanceNoNorm(double* v1, double* v2, int n) {
    int i;
    double nor=0.0;
    for(i=0;i<n;i++) {
        nor+=v1[i]*v2[i];
    }
    nor=1.0-(nor);

    return nor;

}


double** calculateZ(double** X1, int n, int p) {
    int i,j;

    int* group = (int *)malloc(sizeof(int *)*n);
    for(i=0;i<n;i++) group[i]=i;

    double** Z = (double**)malloc(sizeof(double *)*(n-1));
    for(i=0;i<n-1;i++) {
        Z[i]=(double *)malloc(sizeof(double)*3);
        for(j=0;j<3;j++) Z[i][j]=-1.0;
    }

  return Z;
}


double* calculateRelations(double** X1) {
    int i,j;
    double* R = (double *)malloc(sizeof(double)*((int)((graphN)*(graphN-1)/2.0)));
    int cont=0;
    double dist=0.0;

    for(i=0;i<graphN;i++) {
        for(j=i+1;j<graphN;j++) {
            dist = 0.0;
            dist = distance(X1[i],X1[j],P);
            R[cont]=dist;
            if ((R[cont]<0.0) || (R[cont]>1.0)) {
                printf("distance is not correct");
                exit(-1);
            }
            cont++;
        }
    }

    return R;
}


int homegen(int *vect, int n) {

    int max;
    int i;
    int j;
    int jj;

    max=-1;
    for(i=0;i<n;i++) if (vect[i]>max) max=vect[i];
    max=max+1;

    int* mem = (int *)malloc(sizeof(int)*max);
    for(i=0;i<max;i++) mem[i]=-1;

    jj=0;
    for(j=0;j<n;j++) {
        if (mem[vect[j]]==-1) {
            mem[vect[j]]=jj;
            jj++;
        }
    }

    for(j=0;j<n;j++) {
        vect[j]=mem[vect[j]];
    }

    max=-1;
    for(i=0;i<n;i++) if (vect[i]>max) max=vect[i];
    max=max+1;

    return max;
}


int mult_thread(int n, int from, int to, int* seeds, double* conninv, double* sum_weights, double* max_weights, double* val, int* tmp_class) {

    int* cache = (int *)malloc(sizeof(int)*n);
    int* list = (int *)malloc(sizeof(int)*n);
    double* col = (double *)calloc(n,sizeof(double));
    double* col2 = (double *)calloc(n,sizeof(double));
    double* tmpcol;

    int listp;
    double vv;

    int iter;
    int done;
    int i, j, k, k2;

    for(i=from;i<to;i++) {

        iter=0;
        done=0;

        for(j=0;j<n;j++) {
            col[j]=0.0;
            cache[j]=0;
        }
        col[seeds[i]]=1.0;
        cache[seeds[i]]=1;
        list[0]=seeds[i];
        listp=1;

        iter=0;
        while(!done) {

            for(j=0;j<n;j++) col2[j]=0.0;

            for(k=0;k<listp;k++) {
                j=list[k];
                vv = col[j]*(max_weights[j]/sum_weights[j]);
                col2[j] += vv;
                for(k2=0;k2<conn[j];k2++) {
                    if (cache[graph[j][k2]]==0)  {
                        list[listp]=graph[j][k2];
                        listp++;
                        cache[graph[j][k2]]=cache[j]+1;
                    }

                    //vv = col[j]*(conninv_weight[graph[j][k2]]/(double)conn[graph[j][k2]]);
                    vv = col[j]*(graph_weight[j][k2]/sum_weights[j]);
                    col2[graph[j][k2]] += vv;
                }
            }

            tmpcol=col;
            col=col2;
            col2=tmpcol;

            iter++;
            if (iter==3) done=1;
        }


        for(k=0;k<listp;k++) {
            j=list[k];
            //printf("%d %d %f\n",i,j,col[j]);
            if (col[j]>val[j]) {
                val[j]=col[j];
                tmp_class[j]=i;
            }
        }

    }

    //for(i=0;i<n;i++) {
    //    printf("%d %d %f\n",i,tmp_class[i],val[i]);
    //}
    //printf("///////\n");

    // ---

    free(cache);
    free(list);
    free(col);
    free(col2);

    return 0;

}

struct mult_thread_params {
    int n;
    int from;
    int to;
    int* seeds;
    double* conninv;
    double* sum_weights;
    double* max_weights;
    double* val;
    int* tmp_class;
};

void* mult_thread_wrapper(void *data) {

    struct mult_thread_params* p = (struct mult_thread_params *)data;
    mult_thread(p->n, p->from, p->to, p->seeds, p->conninv, p->sum_weights, p->max_weights, p->val, p->tmp_class);
    return 0;
}



int* mult_all_3(int n, int* seeds, int* pp) {

    int* class = (int *)malloc(sizeof(int)*n);
    int i,j;
    int p = (*pp);

    double* conninv = (double *)malloc(sizeof(double)*n);
    double* sum_weights = (double *)malloc(sizeof(double)*n);
    double* max_weights = (double *)malloc(sizeof(double)*n);

    double** val;
    int** tmp_class;

    double sumw;
    double maxw;
    double tmp_val;
    int step;
    int base;

    for(i=0;i<n;i++) {
        class[i]=-1;
        sumw = 0.0;
        maxw = 0.0;
        for(j=0;j<conn[i];j++) {
            sumw += graph_weight[i][j];
            if (maxw < graph_weight[i][j]) {
                maxw = graph_weight[i][j];
            }
        }
        conninv[i] = 1.0/(1.0+(double)conn[i]);
        sum_weights[i] = sumw + maxw;
        max_weights[i] = maxw;
    }

    int num_threads = mNumThreads;

    val = (double**)malloc(sizeof(double*)*num_threads);
    tmp_class = (int**)malloc(sizeof(int*)*num_threads);
    for(i=0;i<num_threads;i++) {
        val[i] = (double *)calloc(n,sizeof(double));
        tmp_class[i] = (int *)calloc(n,sizeof(int));
        for(j=0;j<n;j++) {
            val[i][j] = 0.0;
            tmp_class[i][j] = -1;
        }
    }

    pthread_t* thread_pool = (pthread_t*)malloc(sizeof(pthread_t)*num_threads);
    struct mult_thread_params** thread_pool_params = (struct mult_thread_params**)malloc(sizeof(struct mult_thread_params*)*num_threads);
    struct mult_thread_params* tp;

    step = p / num_threads;
    base = 0;

    for(i=0;i<num_threads;i++) {

        thread_pool_params[i] = (struct mult_thread_params*)malloc(sizeof(struct mult_thread_params));
        tp = thread_pool_params[i];

        tp->n = n;
        tp->from = base;
        if (i==num_threads-1) tp->to = p;
        else tp->to = base + step;
        tp->seeds = seeds;
        tp->conninv = conninv;
        tp->sum_weights = sum_weights;
        tp->max_weights = max_weights;
        tp->val = val[i];
        tp->tmp_class = tmp_class[i];

        if (pthread_create(&thread_pool[i], NULL, mult_thread_wrapper, tp)) {
            printf("Error creating thread!\n");
            exit(-1);
        }

        base += step;
        //the case of no-threads
        //mult_thread(n,0,p,seeds,conninv,sum_weights,max_weights,val[i],tmp_class[i]);
    }

    for(i=0;i<num_threads;i++) {
        pthread_join(thread_pool[i], NULL);
    }

    printf("All threads finished!\n");

    for(i=0;i<n;i++) {
        tmp_val = 0.0;
        for(j=0;j<num_threads;j++) {
            if (val[j][i] > tmp_val) {
                tmp_val = val[j][i];
                class[i] = tmp_class[j][i];
            }
        }
    }


    for(i=0;i<num_threads;i++) {
        free(tmp_class[i]);
        free(val[i]);
        free(thread_pool_params[i]);
    }
    free(tmp_class);
    free(val);


    for(i=0;i<n;i++) {
        if (class[i]==-1) {
            class[i]=(*pp);
            (*pp)=(*pp)+1;
        }
    }

    free(conninv);
    free(sum_weights);
    free(max_weights);

    return class;

}

int min(int x, int y) {
    if (x<y) return x;
    else return y;
}

int max(int x, int y) {
    if (x<y) return y;
    else return x;
}


int build(int* class, int n, int p) {

    int* classtemp=NULL;
    int* classmax=NULL;
    int* classmax_temp=NULL;
    double qqmax=-1.0;
    int maxp=-1;
    int c1;
    int c2;
    int i, j;
    int numlinks=0;


    FILE* fdclass = fopen("./results.class","w");
    FILE* fdmod = fopen("./results.mod","w");
    FILE* fdden = fopen("./results.dendro","w");
    FILE* fdjoin = fopen("./results.join","w");

    FILE* fd_init_class = fopen("./results.init.class","w");


    for(i=0;i<n;i++) {
        if (class[i]>maxp) maxp=class[i];
    }

    maxp = homegen(class,n);

    for(i=0;i<n;i++) {
        if (class[i]<0) {
            printf("error in class: %d %d\n",i,class[i]);
            exit(-1);
        }
    }

    classmax_temp = (int *)malloc(sizeof(int)*n);
    classtemp = (int *)malloc(sizeof(int)*n);
    classmax = (int *)malloc(sizeof(int)*n);
    for(i=0;i<n;i++) classtemp[i]=class[i];

    int* outlinks = (int *)malloc(sizeof(int)*maxp);
    int* inlinks = (int *)malloc(sizeof(int)*maxp);

    for(i=0;i<maxp;i++) {
        outlinks[i]=0;
        inlinks[i]=0;
    }

    for(i=0;i<n;i++) {
        c1=class[i];
        for(j=0;j<conn[i];j++) {
            c2=class[graph[i][j]];

            if (c1==c2) inlinks[c1]++;

            outlinks[c1]++;

            numlinks++;
        }
    }

    double qq = 0.0;
    double* vm = (double *)malloc(sizeof(double)*maxp);

    for(i=0;i<maxp;i++) {
        vm[i] = (((double)inlinks[i]/(double)numlinks) - pow((double)outlinks[i]/(double)numlinks,2.0));
        qq=qq+vm[i];
    }


    fprintf(fdmod,"%f\n",qq);

    // ---------------------------------------------------
    // now the real thing

    double minqq;
    double maxqq;
    int minqqpos;
    int maxqqpos;
    int pos1, pos2;
    double tempqq;
    double newvm=0.0;
    double currnewvm=0.0;
    double err;
    int jj;
    int c3;

    double** Z = (double **)malloc(sizeof(double *)*(graphN-1));
    for(i=0;i<graphN-1;i++) Z[i]=(double *)malloc(sizeof(double)*3);
    for(i=0;i<graphN-1;i++) for(j=0;j<3;j++) Z[i][j]=0.0;
    double* labZ = (double *)malloc(sizeof(double)*maxp);
    for(i=0;i<maxp;i++) labZ[i]=0.0;

    int* tt = (int *)malloc(sizeof(int)*n);
    for(j=0;j<n;j++) tt[j]=0;

    int contz=0;
    int conttt;
    double contzz=0.1;
    //printf("dendrogram\n");
    for(i=0;i<maxp;i++) {

        conttt=0;
        for(j=0;j<n;j++) {
            if (classtemp[j]==i) {
                tt[conttt]=j;
                conttt++;
            }
        }

        if (conttt>1) {
            Z[contz][0]=tt[0];
            Z[contz][1]=tt[1];
            Z[contz][2]=contzz;
            contz=contz+1;
            for(j=2;j<conttt;j++) {
                Z[contz][0]=tt[j];
                Z[contz][1]=contz-1+graphN;
                Z[contz][2]=contzz;
                contz=contz+1;
            }
            labZ[i]=contz-1+graphN;
    }
    else labZ[i]=tt[0];

    //printf("Z %d %d %f\n",(int)Z[contz-1][0],(int)Z[contz-1][1],Z[contz-1][2]);

    }

    contzz=contzz+0.05;
    free(tt);

    int* arv = (int *)malloc(sizeof(int)*maxp);
    int* list = (int *)malloc(sizeof(int)*n);
    int listp;

    pos1=-1;
    for(j=0;j<n;j++) {
        c2 = classtemp[j];
        if (c2!=pos1) {
            for(jj=0;jj<conn[j];jj++) {
                c3 = classtemp[graph[j][jj]];
                if (c2==c3) arv[c2]++;
                if (c3==pos1) arv[c2]=arv[c2]+2;
            }
        }
        else {
            for(jj=0;jj<conn[j];jj++) {
                c3 = classtemp[graph[j][jj]];
                if (c2==c3) arv[c2]++;
            }
        }
    }

    int joinscont=0;
    int** joins = (int **)malloc(sizeof(int*)*maxp);
    for(i=0;i<maxp;i++) {
        joins[i] = (int *)malloc(sizeof(int)*2);
        joins[i][0]=0;
        joins[i][1]=0;
    }

    double modmod = 0.0;

    if (mSave==0) {
        int kk = homegen(classtemp,n);
        fprintf(fdclass,"----\n");
        for(i=0;i<n;i++) {
            fprintf(fdclass,"%d\n",classtemp[i]);
            fprintf(fd_init_class,"%d\n",classtemp[i]);
        }
        modmod = modularity(classtemp,graphN);
        fprintf(fdclass,"%f\n%f\n%d\n",qqmax,modmod,kk);

        printf("modularity: %f\n",modmod);
    }

    //double prev_qqmax = modmod;

    for(i=1;i<maxp;i++) {
        minqq=99999.99;
        maxqq=-999999.99;
        minqqpos=-1;
        maxqqpos=-1;

        //printf("build: %d out of %d\n",i,p);

        for(j=0;j<maxp;j++) {
            if ((outlinks[j]>0) && (vm[j]<minqq)) {
                minqq=vm[j];
                minqqpos=j;
            }
            arv[j]=0;
        }

        pos1 = minqqpos;

        listp=0;
        for(j=0;j<n;j++) {
            if (classtemp[j]==pos1) {
                list[listp]=j;
                listp++;
            }
        }

        for(j=0;j<listp;j++) {
                for(jj=0;jj<conn[list[j]];jj++) {
                    c2 = classtemp[graph[list[j]][jj]];
                    if (c2!=pos1) {
                        arv[c2]+=2;

                    }
                }

        }

        for(j=0;j<maxp;j++) {
            if ((outlinks[j]>0) && (j!=pos1)) {

                tempqq = qq + ((double)arv[j]/(double)numlinks) - (2.0*(((double)outlinks[pos1]/(double)numlinks)*((double)outlinks[j]/(double)numlinks)));

                if (tempqq>maxqq) {
                    maxqq=tempqq;
                    maxqqpos=j;
                    currnewvm=newvm;
                }
                else if (tempqq==maxqq) {
                    if (newvm<currnewvm) {
                        maxqq=tempqq;
                        maxqqpos=j;
                        currnewvm=newvm;
                    }
                }
            }
        }

        pos2=maxqqpos;

        int totpos1=0;
        int totpos2=0;
        for(j=0;j<n;j++) {
            if (classtemp[j]==pos1) totpos1++;
            if (classtemp[j]==pos2) totpos2++;
        }

        joins[joinscont][0]=totpos1;
        joins[joinscont][1]=totpos2;

        joinscont++;

        for(j=0;j<n;j++) {
            if (classtemp[j]==pos2) {
                classtemp[j]=pos1;
            }
        }

        Z[contz][0]=labZ[pos2];
        Z[contz][1]=labZ[pos1];
        Z[contz][2]=contzz;
        labZ[pos1]=contz+graphN;
        labZ[pos2]=contz+graphN;
        contz=contz+1;
        contzz=contzz+0.05;

        tempqq = qq - (vm[pos1]+vm[pos2]);

        err = outlinks[pos1] + outlinks[pos2];

        if (outlinks)

        newvm = ((double)(arv[pos2]+inlinks[pos1]+inlinks[pos2])/(double)numlinks) - pow((double)err/(double)numlinks,2.0);

        vm[pos1]=newvm;
        vm[pos2]=0.0;

        qq = tempqq + vm[pos1];
        fprintf(fdmod,"%f\n",qq);

        if (mSave) {
            for(j=0;j<n;j++) fprintf(fdclass,"%d ",classtemp[j]);
            fprintf(fdclass,"\n");
        }
        else {
            if (qq>qqmax) {
                qqmax=qq;
                for(j=0;j<n;j++) classmax[j]=classtemp[j];
            }
        }


        if (partitions_to_save>0) {

            int curr_num_clusters = (maxp - i);
            int ci;
            for(ci=0;ci<partitions_to_save;ci++) {
                if (partitions_to_save_v[ci]==curr_num_clusters) {
                    char temp_s[255];
                    sprintf(temp_s, "./partition_%d_clusters.class",curr_num_clusters);
                    FILE* fdtemp = fopen(temp_s,"w");

                    int* classtemp_homegen = (int*)malloc(sizeof(int)*n);
                    for(j=0;j<n;j++) classtemp_homegen[j] = classtemp[j];
                    homegen(classtemp_homegen,n);
                    for(j=0;j<n;j++) fprintf(fdtemp,"%d\n",classtemp_homegen[j]);
                    free(classtemp_homegen);
                    fclose(fdtemp);
                }
            }
        }


        inlinks[pos1]=inlinks[pos1]+inlinks[pos2]+arv[pos2];
        inlinks[pos2]=0;
        outlinks[pos1]=outlinks[pos1]+outlinks[pos2];
        outlinks[pos2]=0;

    }

    free(arv);
    free(list);

    // ---------------------------------------------------

    if (mSave==0) {
        int kk = homegen(classmax,n);
        fprintf(fdclass,"----\n");
        for(i=0;i<n;i++) {
            fprintf(fdclass,"%d\n",classmax[i]);
        }
        double modmod = modularity(classmax,graphN);
        fprintf(fdclass,"%f\n%f\n%d\n",qqmax,modmod,kk);

        printf("modularity: %f\n",modmod);
    }

    for(i=0;i<graphN-1;i++) {
        fprintf(fdden,"%d %d %f\n",(int)(Z[i][0]+1),(int)(Z[i][1]+1),Z[i][2]);
    }

    for(i=0;i<maxp-1;i++) free(Z[i]);

    for(i=0;i<maxp-1;i++) {
        fprintf(fdjoin,"%d %d\n",joins[i][0],joins[i][1]);
    }


    free(inlinks);
    free(outlinks);
    free(class);
    free(classtemp);
    free(vm);
    free(labZ);
    free(Z);
    free(classmax);

    fclose(fdclass);
    fclose(fdmod);
    fclose(fdden);
    fclose(fdjoin);
    fclose(fd_init_class);


  return 0;

}



int main(int argc, char* argv[]) {

    int i;

    srand(-((time(NULL)%100)+1));

    long t1 = time(NULL);
    struct timeb tt1, tt2;
    ftime(&tt1);
    int* class;


    if (argc!=6) {
        printf("\nUsage:\n\n./clusterPBD filename type distance msave num_cores\n\n");
        printf("   filename (network in pajek format)\n");
        printf("   type {0,1} (pajek format: 0, pajek compact format: 1)\n");
        printf("   distance (for the initial seed)\n");
        printf("   msave {0,1} (don't save intermediate results: 0, save them: 1)\n");
        printf("   num_cores number of threads, one per core recommended\n");
        printf("\n");
        printf("Examples:\n\n");
        printf("   ./clusterPBD data/zachary.net 0 2 0 4\n");
        printf("   ./clusterPBD data/erdos02b.net 1 2 0 4\n");
        printf("\n");
        exit(-1);
    }

    int type = atoi(argv[2]);
    D = atoi(argv[3]);

    mSave = atoi(argv[4]);
    mNumThreads = atoi(argv[5]);


    if (type==0) load(argv[1]);
    else loadShort(argv[1]);

    printf("Network loaded successfully\n");

    for(i=0;i<graphN;i++) {
        if (conn[i]==0) {
            printf("The network is not fully connected (problem in node %d)\n",i);
            exit(-1);
        }
    }

    // part of the cliqz hack
    partitions_to_save = 8;
    if (partitions_to_save>0) {
        partitions_to_save_v = (int*)malloc(sizeof(int)*partitions_to_save);
        partitions_to_save_v[0] = 9;
        partitions_to_save_v[1] = 50;
        partitions_to_save_v[2] = 100;
        partitions_to_save_v[3] = 200;
        partitions_to_save_v[4] = 500;
        partitions_to_save_v[5] = 1000;
        partitions_to_save_v[6] = 5000;
        partitions_to_save_v[7] = 10000;
    }


    // true if normal behavior, otherwise the seed stage is loaded from file
    if (1==1) {
        int* seeds = matrixSeedMaxConn(graphN,&P);
        printf("Seed calculated %d\n",P);

        class = mult_all_3(graphN,seeds,&P);
        printf("Spectral phase done %d\n",P);

    }
    else {
        printf("You are using the hacked version only for Wreduced.pajek.net!\n");
        class = load_initial_partition("./Wreduced_initial_partition_after_spectral.part", &P);
        printf("loaded!!! %d %d\n",P,graphN);
    }

    build(class,graphN,P);

    long t2 = time(NULL);
    ftime(&tt2);

    printf("execution time: %ld %ld %lds\n",t2,t1,t2-t1);
    printf("N: %d E: %d\n",graphN,graphE);
    return 0;
}
