#include <stdio.h>                  /*  Marr-Hildreth.c  (or marrh.c) */
#include <math.h>
#define  PICSIZE 256
#define  MAXMASK 100
#define  TRUE 1
#define  FALSE 0

         int      in_pic[PICSIZE][PICSIZE];
         double   x_mask[MAXMASK][MAXMASK];
         double   y_mask[MAXMASK][MAXMASK];
         double   x_conv[PICSIZE][PICSIZE];
         double   y_conv[PICSIZE][PICSIZE];
         double   magnitude[PICSIZE][PICSIZE];
         double   cand[PICSIZE][PICSIZE];
         double   final[PICSIZE][PICSIZE];
         double   histogram[PICSIZE];
         int      slope;
         int      rows, cols;

         // function declarations
         void readHeader(FILE* imgFin);      // function to read .pgm header
         int isComment(char *line);          // for readHeader() function


// main
main(argc,argv)
int argc;
char **argv;
{
   int      i, j, p, q, s, t, mr, centx, centy;
   double   x_maskval, y_maskval, sum_x, sum_y, sig, maxival, minival, maxval, ZEROTOL, percent;
   FILE     *fo1, *fo2, *fo3, *fp1, *fopen();
   char     *foobar;

   // IN file
   argc--; argv++;
   foobar = *argv;
   fp1 = fopen("garb34.pgm","rb");
   readHeader(fp1);

   // OUT file 1 - magnitude
   argc--; argv++;
   foobar = *argv;
   fo1=fopen("canny_mag.pgm","wb");

   fprintf(fo1, "P5\n");
   fprintf(fo1, "%d %d\n", rows, cols);
   fprintf(fo1, "255\n");

   // OUT file 2 - peaks
   argc--; argv++;
   foobar = *argv;
   fo2=fopen("canny_peaks.pgm","wb");

   fprintf(fo2, "P5\n");
   fprintf(fo2, "%d %d\n", rows, cols);
   fprintf(fo2, "255\n");

   // OUT file 3 - final
   argc--; argv++;
   foobar = *argv;
   fo3=fopen("canny_final.pgm","wb");

   fprintf(fo3, "P5\n");
   fprintf(fo3, "%d %d\n", rows, cols);
   fprintf(fo3, "255\n");

   // sig = atof(foobar);
   sig = 1;       

   // argc--; argv++;
   // foobar = *argv;
   // ZEROTOL = atof(foobar);

   // flexible size masks
   mr = (int)(sig * 3);
   centx = (MAXMASK / 2);
   centy = (MAXMASK / 2);

   for (i=0;i<256;i++)
   { for (j=0;j<256;j++)
      {
         in_pic[i][j] =  getc (fp1);
         in_pic[i][j] &= 0377;
      }
   }

   // first derivitives
   for (p=-mr;p<=mr;p++)
   {  for (q=-mr;q<=mr;q++)
      {
         x_maskval = (-1 * q ) * (exp(-1 * (((p * p) + (q * q)) / (2 * sig * sig))));
         x_mask[p+centy][q+centx] = x_maskval;
         
         y_maskval = (-1 * p ) * (exp(-1*(((p*p)+(q*q))/(2*(sig*sig)))));
         y_mask[p+centy][q+centx] = y_maskval;
      }
   }


   //convolution
   for (i=0;i<=255;i++)
   { for (j=0;j<=255;j++)
      {
         sum_x = 0;
         sum_y = 0;

         for (p=-mr;p<=mr;p++)
         {
            for (q=-mr;q<=mr;q++)
            {
               sum_x += in_pic[i+p][j+q] * x_mask[p+centy][q+centx];
               sum_y += in_pic[i+p][j+q] * y_mask[p+centy][q+centx];
            }
         }
         x_conv[i][j] = sum_x;
         y_conv[i][j] = sum_y;
      }
   }

   // from sobel.c
   // Part 1
   maxval  = 0;
   for (i=mr;i<256-mr;i++)
   { for (j=mr;j<256-mr;j++)
      {
         magnitude[i][j]=sqrt((double)((x_conv[i][j]*x_conv[i][j]) + (y_conv[i][j]*y_conv[i][j])));
         if (magnitude[i][j] > maxival)
            maxival = magnitude[i][j];
      }
   }


   for (i=0;i<256;i++)
   { for (j=0;j<256;j++)
      {
         magnitude[i][j] = (magnitude[i][j] / maxival) * 255;
         fprintf(fo1,"%c",(char)((int)(magnitude[i][j])));

      }
   }

   // peaks - from module 2 part 2 lecture slide #24 (../reulectureOnline1.pptx)
   // Part 2
   for(i=mr;i<256-mr;i++){
      for(j=mr;j<256-mr;j++){
         if((x_conv[i][j]) == 0.0) {
            x_conv[i][j] = .00001;
         }
         slope = y_conv[i][j]/x_conv[i][j];
         if( (slope <= .4142)&&(slope > -.4142)){
            if((magnitude[i][j] > magnitude[i][j-1])&&(magnitude[i][j] > magnitude[i][j+1])){
               cand[i][j] = 255;
            }
         }
         else if( (slope <= 2.4142)&&(slope > .4142)){
            if((magnitude[i][j] > magnitude[i-1][j-1])&&(magnitude[i][j] > magnitude[i+1][j+1])){
               cand[i][j] = 255;
            }
         }
         else if( (slope <= -.4142)&&(slope > -2.4142)){
            if((magnitude[i][j] > magnitude[i+1][j-1])&&(magnitude[i][j] > magnitude[i-1][j+1])){
               cand[i][j] = 255;
            }
         }else{
            if((magnitude[i][j] > magnitude[i-1][j])&&(magnitude[i][j] > magnitude[i+1][j])){
               cand[i][j] = 255;
            }
         }
     }
  }


   for (i=0;i<256;i++)
   { for (j=0;j<256;j++)
      {
         fprintf(fo2,"%c",(char)((int)(cand[i][j])));
      }
   }

   // Part 4
   printf("Enter a percentage (0.5%% as 0.005): ");
   scanf("%lf", &percent);

   int cutoff = percent * rows * cols;
   int histogramSize;
   int HI, LO, areaOfTops;

   // histogram of scaled magnitudes
   for(i=mr;i<256-mr;i++)
      for(j=mr;j<256-mr;j++)
         histogram[(int)magnitude[i][j]]++;
   histogramSize = sizeof(histogram) / sizeof(histogram[0]);

   for(i=histogramSize;i>=0;i--){
      areaOfTops += histogram[i];
      if(areaOfTops > cutoff){
         HI = i;
         break;
      }
   }
   LO = .35*HI;

   printf("rows=%d\tcols=%d\tpercent=%d\n",rows,cols,percent);
   printf("cutoff=%d\tareaOfTops=%d\n",cutoff,areaOfTops);
   printf("HI=%d\tLO=%d\n", HI,LO);
   printf("histogramSize=%d\n",histogramSize);


   // Part 3
   for(i=mr; i<256-mr;i++)
   { for(j=mr;j<256-mr;j++)
      {
         if (cand[i][j] == 255){
            if (magnitude[i][j] > HI){
                  cand[i][j] = 0;
                  final[i][j] = 255;
            }
            else if (magnitude[i][j]< LO){
               cand[i][j] = 0;
               final[i][j] = 0;
            }
         }
      }
   }

   int moretodo = TRUE;

   while(moretodo){
      moretodo = FALSE;
      for(i=mr;i<256-mr;i++)
      { for(j=mr;j<256-mr;j++)
         {
            if(cand[i][j] == 255){
               for(p=-1;p<=1;p++)
               { for(q=-1;q<=1;q++)
                  {
                     if(final[i+p][j+q] == 255){
                        cand[i][j] = 0;
                        final[i][j] = 255;
                        moretodo = TRUE;
                     }
                  }

               }
            }
         }
      }
   }

   for (i=0;i<256;i++)
   { for (j=0;j<256;j++)
      {
         fprintf(fo3,"%c",(char)((int)(final[i][j])));
      }

   }
} // end main

// function to read .pgm header
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

// for readHeader() function
int isComment(char *line)
{
        if(line[0] == '#')
        return(TRUE);

        return(FALSE);
}