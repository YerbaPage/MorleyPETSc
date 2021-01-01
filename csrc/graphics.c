/*
 *		graphics.c 
 *
 *		Matrix-Solver Community Project
 *
 *------------------------------------------------------
 *
 *		Created by Chensong Zhang on 03/29/2009.
 *		Copyright 2009 PSU. All rights reserved. 
 *
 *------------------------------------------------------
 *
 */

/*! \file graphics.c
 *  \brief Functions for graphical output
 */

#include <stdio.h>
#include <math.h>

/**
 *  NAME
 *
 *  write_bmp16 - write 16-color raster image in BMP file format
 *
 *  DESCRIPTION
 *
 *  The routine write_bmp16 writes 16-color raster image in
 *  uncompressed BMP file format (Windows bitmap) to a binary file whose
 *  name is specified by the character string fname.
 *
 *  The parameters m and n specify, respectively, the number of rows and
 *  the numbers of columns (i.e. height and width) of the raster image.
 *
 *  The character array map has m*n elements. Elements map[0, ..., n-1]
 *  correspond to the first (top) scanline, elements map[n, ..., 2*n-1]
 *  correspond to the second scanline, etc.
 *
 *  Each element of the array map specifies a color of the corresponding
 *  pixel as 8-bit binary number XXXXIRGB, where four high-order bits (X)
 *  are ignored, I is high intensity bit, R is red color bit, G is green
 *  color bit, and B is blue color bit. Thus, all 16 possible colors are
 *  coded as following hexadecimal numbers:
 *
 *     0x00 = black         0x08 = dark gray
 *     0x01 = blue          0x09 = bright blue
 *     0x02 = green         0x0A = bright green
 *     0x03 = cyan          0x0B = bright cyan
 *     0x04 = red           0x0C = bright red
 *     0x05 = magenta       0x0D = bright magenta
 *     0x06 = brown         0x0E = yellow
 *     0x07 = light gray    0x0F = white
 *
 *  RETURNS
 *
 *  If no error occured, the routine returns zero; otherwise, it prints
 *  an appropriate error message and returns non-zero. 
 *  This code is modified from graphical output in 
 *         GLPK (GNU Linear Programming Kit).
 *
 *  Modified by Chensong Zhang
 *
 *  Note: 
 *
 *  GLPK is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  GLPK is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
 */

static void put_byte(FILE *fp, int c)
{ fputc(c, fp);
	return;
}

static void put_word(FILE *fp, int w)
{ /* big endian */
	put_byte(fp, w);
	put_byte(fp, w >> 8);
	return;
}

static void put_dword(FILE *fp, int d)
{ /* big endian */
	put_word(fp, d);
	put_word(fp, d >> 16);
	return;
}

/*! 
 * \fn int write_bmp16(const char *fname, int m, int n, const char map[])
 * \brief Write to BMP file
 * \param fname char for filename 
 * \param m number of pixels for height
 * \param n number of pixels for weight
 * \param map picture for BMP  
 * \return 1 if succeed, 0 if fail
 */
int write_bmp16(const char *fname, int m, int n, const char map[])
{     FILE *fp;
	int offset, bmsize, i, j, b, ret = 1;
	if (!(1 <= m && m <= 32767))
		printf("write_bmp16: m = %d; invalid height\n", m);
	if (!(1 <= n && n <= 32767))
		printf("write_bmp16: n = %d; invalid width\n", n);
	fp = fopen(fname, "wb");
	if (fp == NULL)
	{  printf("write_bmp16: unable to create `%s'\n", fname);
		ret = 0;
		goto fini;
	}
	offset = 14 + 40 + 16 * 4;
	bmsize = (4 * n + 31) / 32;
	/* struct BMPFILEHEADER (14 bytes) */
	/* UINT bfType */          put_byte(fp, 'B'), put_byte(fp, 'M');
	/* DWORD bfSize */         put_dword(fp, offset + bmsize * 4);
	/* UINT bfReserved1 */     put_word(fp, 0);
	/* UNIT bfReserved2 */     put_word(fp, 0);
	/* DWORD bfOffBits */      put_dword(fp, offset);
	/* struct BMPINFOHEADER (40 bytes) */
	/* DWORD biSize */         put_dword(fp, 40);
	/* LONG biWidth */         put_dword(fp, n);
	/* LONG biHeight */        put_dword(fp, m);
	/* WORD biPlanes */        put_word(fp, 1);
	/* WORD biBitCount */      put_word(fp, 4);
	/* DWORD biCompression */  put_dword(fp, 0 /* BI_RGB */);
	/* DWORD biSizeImage */    put_dword(fp, 0);
	/* LONG biXPelsPerMeter */ put_dword(fp, 2953 /* 75 dpi */);
	/* LONG biYPelsPerMeter */ put_dword(fp, 2953 /* 75 dpi */);
	/* DWORD biClrUsed */      put_dword(fp, 0);
	/* DWORD biClrImportant */ put_dword(fp, 0);
	/* struct RGBQUAD (16 * 4 = 64 bytes) */
	/* CGA-compatible colors: */
	/* 0x00 = black */         put_dword(fp, 0x000000);
	/* 0x01 = blue */          put_dword(fp, 0x000080);
	/* 0x02 = green */         put_dword(fp, 0x008000);
	/* 0x03 = cyan */          put_dword(fp, 0x008080);
	/* 0x04 = red */           put_dword(fp, 0x800000);
	/* 0x05 = magenta */       put_dword(fp, 0x800080);
	/* 0x06 = brown */         put_dword(fp, 0x808000);
	/* 0x07 = light gray */    put_dword(fp, 0xC0C0C0);
	/* 0x08 = dark gray */     put_dword(fp, 0x808080);
	/* 0x09 = bright blue */   put_dword(fp, 0x0000FF);
	/* 0x0A = bright green */  put_dword(fp, 0x00FF00);
	/* 0x0B = bright cyan */   put_dword(fp, 0x00FFFF);
	/* 0x0C = bright red */    put_dword(fp, 0xFF0000);
	/* 0x0D = bright magenta */ put_dword(fp, 0xFF00FF);
	/* 0x0E = yellow */        put_dword(fp, 0xFFFF00);
	/* 0x0F = white */         put_dword(fp, 0xFFFFFF);
	/* pixel data bits */
	b = 0;
	for (i = m - 1; i >= 0; i--) {  
		for (j = 0; j < ((n + 7) / 8) * 8; j++) {  
			b <<= 4;
			b |= (j < n ? map[i * n + j] & 15 : 0);
			if (j & 1) put_byte(fp, b);
		}
	}
	fflush(fp);
	if (ferror(fp)) {  
		printf("write_bmp16: write error on `%s'\n",fname);
		ret = 0;
	}
fini: if (fp != NULL) fclose(fp);
	return ret;
}


/**
 * \struct grid2d
 * \brief 2d grid data structure
 *
 * The grid2d structure is simply a list of triangles, edges and vertices. 
 *			edge i has 2 vertices e[i], 
 *			triangle i has 3 edges s[i], 3 vertices t[i]
 *			vertex i has two coordinates p[i]
 */
typedef struct grid2d {
	double (*p)[2];	/**< Coordinates of vertices */
	int (*e)[2];			/**< Vertices of edges */
	int (*t)[3];			/**< Vertices of triangles */
	int (*s)[3];			/**< Edges of triangles */
	int *pdiri;			/**< Boundary flags (0 <=> interior point) */
	int *ediri;			/**< Boundary flags (0 <=> interior edge) */
	
	int *pfather;		/**< Father point or edge */
	int *efather;		/**< Father edge or triangle */
	int *tfather;		/**< Father triangle */
	
	int vertices;		/**< Number of grid points */
	int edges;				/**< Number of edges */
	int triangles;		/**< Number of triangles */
} grid2d;
typedef grid2d *pgrid2d;
typedef const grid2d *pcgrid2d;

/*! 
 * \fn void output_grid2d(pgrid2d pg, int level)
 * \brief Output grid to a EPS file
 * \param pg pointer to grid in 2d
 * \param level number of level  
 */
void output_grid2d(pgrid2d pg, int level)
{
	FILE *datei;
	char buf[120];
	int i;
	double xmid,ymid,xc,yc;
	
	sprintf(buf,"Grid_ref_level%d.eps",level);
	datei = fopen(buf,"w");
	fprintf(datei, "%%!PS-Adobe-2.0-2.0 EPSF-2.0\n");
	fprintf(datei, "%%%%BoundingBox: 0 0 550 550\n");
	fprintf(datei, "25 dup translate\n");
	fprintf(datei, "%f setlinewidth\n",0.2);
	fprintf(datei, "/Helvetica findfont %f scalefont setfont\n",64.0*pow(0.5,level));
	fprintf(datei, "/b{0 setgray} def\n");
	fprintf(datei, "/r{1.0 0.6 0.6  setrgbcolor} def\n");
	fprintf(datei, "/u{0.1 0.7 0.1  setrgbcolor} def\n");
	fprintf(datei, "/d{0.1 0.1 1.0  setrgbcolor} def\n");
	fprintf(datei, "/cs{closepath stroke} def\n");
	fprintf(datei, "/m{moveto} def\n");
	fprintf(datei, "/l{lineto} def\n");
	
	fprintf(datei,"b\n");
	for(i=0; i<pg->triangles; i++){
		xc = (pg->p[pg->t[i][0]][0]+pg->p[pg->t[i][1]][0]+pg->p[pg->t[i][2]][0])*150.0;
		yc = (pg->p[pg->t[i][0]][1]+pg->p[pg->t[i][1]][1]+pg->p[pg->t[i][2]][1])*150.0;
		
		xmid = pg->p[pg->t[i][0]][0]*450.0;
		ymid = pg->p[pg->t[i][0]][1]*450.0;
		fprintf(datei,"%.1f %.1f m ",0.9*xmid+0.1*xc,0.9*ymid+0.1*yc);
		xmid = pg->p[pg->t[i][1]][0]*450.0;
		ymid = pg->p[pg->t[i][1]][1]*450.0;
		fprintf(datei,"%.1f %.1f l ",0.9*xmid+0.1*xc,0.9*ymid+0.1*yc);
		xmid = pg->p[pg->t[i][2]][0]*450.0;
		ymid = pg->p[pg->t[i][2]][1]*450.0;
		fprintf(datei,"%.1f %.1f l ",0.9*xmid+0.1*xc,0.9*ymid+0.1*yc);
		fprintf(datei,"cs\n");
	}
	fprintf(datei,"r\n");
	for(i=0; i<pg->vertices; i++){
		xmid = pg->p[i][0]*450.0;
		ymid = pg->p[i][1]*450.0;
		fprintf(datei,"%.1f %.1f m ",xmid,ymid);
		fprintf(datei,"(%d) show\n ",i);
	}
	fprintf(datei,"u\n");
	for(i=0; i<pg->edges; i++){
		xmid = 0.5*(pg->p[pg->e[i][0]][0]+pg->p[pg->e[i][1]][0])*450.0;
		ymid = 0.5*(pg->p[pg->e[i][0]][1]+pg->p[pg->e[i][1]][1])*450.0;
		fprintf(datei,"%.1f %.1f m ",xmid,ymid);
		fprintf(datei,"(%d) show\n ",i);
		
		xmid = pg->p[pg->e[i][0]][0]*450.0;
		ymid = pg->p[pg->e[i][0]][1]*450.0;
		fprintf(datei,"%.1f %.1f m ",xmid,ymid);
		xmid = pg->p[pg->e[i][1]][0]*450.0;
		ymid = pg->p[pg->e[i][1]][1]*450.0;
		fprintf(datei,"%.1f %.1f l ",xmid,ymid);
		fprintf(datei,"cs\n");
	}
	fprintf(datei,"d\n");
	for(i=0; i<pg->triangles; i++){
		xmid = (pg->p[pg->t[i][0]][0]+pg->p[pg->t[i][1]][0]+pg->p[pg->t[i][2]][0])*150.0;
		ymid = (pg->p[pg->t[i][0]][1]+pg->p[pg->t[i][1]][1]+pg->p[pg->t[i][2]][1])*150.0;
		fprintf(datei,"%.1f %.1f m ",xmid,ymid);
		fprintf(datei,"(%d) show\n ",i);
	}
	fprintf(datei, "showpage\n");
	fclose(datei);	
}
