/*
 *  MrBayes 3.1.2
 *
 *  copyright 2002-2005
 *
 *  John P. Huelsenbeck
 *  Section of Ecology, Behavior and Evolution
 *  Division of Biological Sciences
 *  University of California, San Diego
 *  La Jolla, CA 92093-0116
 *
 *  johnh@biomail.ucsd.edu
 *
 *	Fredrik Ronquist
 *  Paul van der Mark
 *  School of Computational Science
 *  Florida State University
 *  Tallahassee, FL 32306-4120
 *
 *  ronquist@scs.fsu.edu
 *  paulvdm@scs.fsu.edu
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details (www.gnu.org).
 *
 */

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

int      CopyResults (FILE *toFile, char *fromFileName, int lastGen);
int      CopyTreeResults (FILE *toFile, char *fromFileName, int lastGen, int *treeNum);
safeLong FirstTree (FILE *fp, char *lineBuf, int longestLine);
int      Flip01 (int x);
void     FlipBits (safeLong *partition, int length, safeLong *mask);
void     FlipOneBit (int n, safeLong *p);
int      GetNameFromNamestring (char *name, const char *nameString, int index);
int      GetNameFromString (char *s, char *tkn, int n);
int      IsBitSet (int i, safeLong *bits);
int      IsPartNested (safeLong *smaller, safeLong *larger, int length);
int      IsPartCompatible (safeLong *smaller, safeLong *larger, int length);
safeLong LastBlock (FILE *fp, char *lineBuf, int longestLine);
int		 LineTermType (FILE *fp);
int      LongestLine (FILE *fp);
void     MrBayesPrint (char *format, ...);
void     MrBayesPrintf (FILE *f, char *format, ...);
FILE    *OpenBinaryFileR (char *name);
FILE 	*OpenTextFileA (char *name);
FILE    *OpenTextFileR (char *name);
FILE 	*OpenTextFileW (char *name);
int      SafeFclose(FILE **fp);
void    *SafeMalloc(size_t s);
char    *SafeStrcat(char **target, const char *source);
void     StripComments (char *s);
