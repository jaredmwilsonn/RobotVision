#include <stdio.h>                          /* Sobel.c */
#include <math.h>
#define TRUE 1
#define FALSE 0

        double max_ival;
        int rows, cols;

        //function declarations
        void threshold(double img_vals[256][256]);
        void also_threshold(double img_vals[256][256]);         // identical to threshold() function but for output2
        void readHeader(FILE* imgFin);                          // function to read .pgm header
        int isComment(char *line);                              // for readHeader() function



void threshold(double img_vals[256][256])
{
        FILE* fo3;
        fo3 = fopen("sobel_out1.pgm", "wb");
        int i, j;

        // .pgm header
        fprintf(fo3,"P5\n");
        fprintf(fo3,"%d %d\n", rows, cols);
        fprintf(fo3,"255\n");


        for (i=0;i<256;i++)
        { 
                for (j=0;j<256;j++)
                {
                        img_vals[i][j] = (img_vals[i][j] / max_ival) * 255;

                        if ((int)(img_vals[i][j]) < 40)
                                img_vals[i][j] = 0;
                        else 
                                img_vals[i][j] = 255;

                        fprintf(fo3,"%c",(char)((int)(img_vals[i][j])));
                }
        }
}

// identical to threshold function but for output2
void also_threshold(double img_vals[256][256])
{

        FILE* fo2;
        fo2 = fopen("sobelout_2.pgm", "wb");
        int i, j;

        // .pgm header
        fprintf(fo2,"P5\n");
        fprintf(fo2,"%d %d\n", rows, cols);
        fprintf(fo2,"255\n");


        for (i=0;i<256;i++)
        { 
                for (j=0;j<256;j++)
                {
                        img_vals[i][j] = (img_vals[i][j] / max_ival) * 255;

                        if ((int)(img_vals[i][j]) < 40)
                                img_vals[i][j] = 0;
                        else 
                                img_vals[i][j] = 255;


                        fprintf(fo2,"%c",(char)((int)(img_vals[i][j])));
                }
        }
}


// read infile && determine # of rows and cols
void readHeader(FILE* imgFin)
{
        int readID = FALSE;
        int readSize = FALSE;
        int readMaxVal = FALSE;
        char line[100];

        while(!(readID && readSize && readMaxVal))
        {
        fgets(line, 100, imgFin);

        if((strlen(line) == 0) || (strlen(line) == 1))
        continue;

        if(isComment(line))
        continue;

        if(!readID)
                readID = TRUE;
        else if(!readSize)
        {
                sscanf(line, "%d %d", &rows, &cols);
                readSize = TRUE;
        }
        else if(!readMaxVal)
                readMaxVal = TRUE;
        }
}

int isComment(char *line)
{
        if(line[0] == '#')
        return(TRUE);

        return(FALSE);
}

main(argc,argv)
int argc;
char **argv;
{
        int i,j,p,q,mr,sum1,sum2;
        FILE  *fp1, *fopen();
        double ival[256][256], ival_1[256][256];
        int pic[256][256];
        int outpicx[256][256];
        int outpicy[256][256];
        int maskx[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
        int masky[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};

        fp1=fopen("garb34.pgm","rb");
        readHeader(fp1);


        for (i=0;i<256;i++)
        { 
                for (j=0;j<256;j++)
                {
                        pic[i][j]  =  getc (fp1);
                        pic[i][j]  &= 0377;
                }
        }

        mr = 1;
        for (i=mr;i<256-mr;i++)
        { 
                for (j=mr;j<256-mr;j++)
                {
                        sum1 = 0;
                        sum2 = 0;
                        for (p=-mr;p<=mr;p++)
                        {
                                for (q=-mr;q<=mr;q++)
                                {
                                        sum1 += pic[i+p][j+q] * maskx[p+mr][q+mr];
                                        sum2 += pic[i+p][j+q] * masky[p+mr][q+mr];
                                }
                        }
                        outpicx[i][j] = sum1;
                        outpicy[i][j] = sum2;
                }
        }

        max_ival = 0;
        for (i=mr;i<256-mr;i++)
        { 
                for (j=mr;j<256-mr;j++)
                {
                        ival[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                                (outpicy[i][j]*outpicy[i][j])));
                        ival_1[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                                (outpicy[i][j]*outpicy[i][j])));
                        if (ival[i][j] > max_ival)
                                max_ival = ival[i][j];

                }
        }

        FILE *fo1;
        fo1=fopen("sobel_mag.pgm","wb");


        //.pgm header
        fprintf(fo1,"P5\n");
        fprintf(fo1,"%d %d\n", 256, 256);
        fprintf(fo1,"255\n");


        threshold(ival_1);
        for (int i=0;i<256;i++)
        { 
                for (int j=0;j<256;j++)
                {
                        ival[i][j] = (ival[i][j] / max_ival) * 255;            
                        fprintf(fo1,"%c",(char)((int)(ival[i][j])));
                }
        }

        also_threshold(ival);

}
