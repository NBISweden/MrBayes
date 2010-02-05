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

/* id-string for ident, do not edit: cvs will update this string */
const char utilsID[]="$Id: utils.c,v 1.7 2009/01/05 14:09:49 ronquist Exp $";

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "mb.h"
#include "globals.h"
#include "utils.h"





/* CopyResults: copy results from one file to another*/
int CopyResults (FILE *toFile, char *fromFileName, int lastGen)
{
    int     longestLine;
    char    *strBuf, *strCpy, *word;
    FILE    *fromFile;

    if ((fromFile = OpenBinaryFileR(fromFileName)) == NULL)
        return ERROR;

    longestLine = LongestLine(fromFile)+10;
    SafeFclose(&fromFile);
    strBuf = (char *) calloc (2*(longestLine+2),sizeof(char));
    strCpy = strBuf + longestLine + 2;

    if ((fromFile = OpenTextFileR(fromFileName)) == NULL)
        return ERROR;
    
    while (fgets(strBuf,longestLine,fromFile)!=NULL)
        {
        strncpy(strCpy,strBuf,longestLine);
        word = strtok(strCpy," ");
        /* atoi returns 0 when word is not integer number */
        if (atoi(word)>lastGen)
            break;
        fprintf (toFile,"%s",strBuf);
        fflush (toFile);
        }
    
    SafeFclose(&fromFile);
    free(strBuf);
    return (NO_ERROR);
}




/* CopyTreeResults: copy tree results from one file to another*/
int CopyTreeResults (FILE *toFile, char *fromFileName, int lastGen, int *numTrees)
{
    int     longestLine;
    char    *strBuf, *strCpy, *word;
    FILE    *fromFile;
    
    (*numTrees) = 0;

    if ((fromFile = OpenBinaryFileR(fromFileName)) == NULL)
        return ERROR;

    longestLine = LongestLine(fromFile)+10;
    SafeFclose(&fromFile);
    strBuf = (char *) calloc (2*(longestLine+2),sizeof(char));
    strCpy = strBuf + longestLine + 2;

    if ((fromFile = OpenTextFileR(fromFileName)) == NULL)
        return ERROR;
    
    while (fgets(strBuf,longestLine,fromFile)!=NULL)
        {
        strncpy(strCpy,strBuf,longestLine);
        word = strtok(strCpy," ");
        if (strcmp(word,"tree")==0)
            {
            word = strtok(NULL," ");
            /* atoi returns 0 when word is not integer number,
               4 is offset to get rid of "rep." in tree name */
            if (atoi(word+4)>lastGen)
                break;
            (*numTrees)++;
            fprintf (toFile,"%s",strBuf);
            }
        else if (*numTrees == 0)   /* do not print the end statement */
            fprintf (toFile,"%s",strBuf);
        fflush (toFile);
        }
        
    SafeFclose(&fromFile);
    free(strBuf);
    return (NO_ERROR);
}





/* FirstTree: Return file position of first tree after current position */
safeLong FirstTree (FILE *fp, char *lineBuf, int longestLine)
{
	safeLong	firstTree;
	char		*word;
	
	do {
		firstTree = ftell(fp);
		if ((fgets (lineBuf, longestLine, fp)) == NULL)
			return 0;
		word = strtok (lineBuf, " ");
		} while (strcmp(word,"tree")!=0);

	return (firstTree);
}





int Flip01 (int x)

{

	if (x == 0)
		return (1);
	else
		return (0);
		
}





void FlipBits (safeLong *partition, int length, safeLong *mask)

{

	int			i;
	
	for (i=0; i<length; i++)
		{
		partition[i] ^= mask[i];
		}
}





/*-----------------------------------------------------------------
|
|	FlipOneBit: flip bit n in safeLong *p
|
------------------------------------------------------------------*/
void FlipOneBit (int n, safeLong *p)

{

	safeLong		x;

	p += n/nBitsInALong;
	x = 1 << (n % nBitsInALong);
	(*p) ^= x;

}





int GetNameFromNamestring (char *name, const char *nameString, int index)

{

	int		i;
	const char   *p;
	
	for (i=0, p=&nameString[0]; (*p)!='\0'; p++)
		{
		if (i == index)
			break;
		if (*p == '|')
			i++;
		}

	if (i < index)
		return (ERROR);

	while (*p!='|' && *p!='\0')
		(*name++) = (*p++);

	*name = '\0';

	return (NO_ERROR);
	
}





int GetNameFromString (char *s, char *tkn, int n)

{

	int		i, j, startI=0, numPrev;
	   	
	i = numPrev = 0;
	while (s[i] != '\0')
		{
		if (s[i] == '|')
			numPrev++;
		i++;
		}
	if (numPrev < n)
		{
		return (ERROR);
		}
		
	if (n == 1)
		startI = 0;
	else
		{
		i = j = 0;
		while (s[i] != '\0')
			{
			if (s[i] == '|')
				j++;
			i++;
			if (j == n - 1)
				{
				startI = i;
				break;
				}
			}
		}
		
	if (s[startI] == '\0')
		{
		MrBayesPrint ("%s   String is too full\n", spacer);
		return (ERROR);
		}
	
	i = startI;
	j = 0;
	while(s[i] != '\0' && s[i] != '|')
		{
		tkn[j++] = s[i++];
		if (s[i] == '\0')
			{
			MrBayesPrint ("%s   String is too full\n", spacer);
			return (ERROR);
			}
		}
	tkn[j] = '\0';

	return (NO_ERROR);
	
}





int IsBitSet (int i, safeLong *bits)

{

	safeLong		x;

	bits += i / nBitsInALong;

	x = 1 << (i % nBitsInALong);

	if ((*bits) & x)
		return (YES);
	else
		return (NO);
		
}





int IsPartCompatible (safeLong *smaller, safeLong *larger, int length)

{

	int i;

	for (i=0; i<length; i++)
		if ((smaller[i]&larger[i]) != 0)
			break;

	if (i != length)	/* potentially incompatible */
		{
		for (i=0; i<length; i++)
			if ((smaller[i]|larger[i]) != larger[i])
				break;
		}
		
	if (i == length)	/* passed either one of the tests */
		return 1;
	else
		return 0;
}





int IsPartNested (safeLong *smaller, safeLong *larger, int length)

{

	int i;

	for (i=0; i<length; i++)
		if ((smaller[i] | larger[i]) != larger[i])
			break;
		
	if (i == length)
		return 1;
	else
		return 0;

}





/* LastBlock: Return file position of last block in file */
safeLong LastBlock (FILE *fp, char *lineBuf, int longestLine)
{
	safeLong	lastBlock;
	char	*word;
	
	lastBlock = 0L;
	rewind (fp);

	while ((fgets (lineBuf, longestLine, fp)) != NULL)
		{
		word = strtok (lineBuf, " ");
		if (strcmp (word, "begin") == 0)
			lastBlock = ftell (fp);
		}

	return lastBlock;
}





int LineTermType (FILE *fp)

{

	int			ch, nextCh, term;

	term = LINETERM_UNIX;	/* default if no line endings are found */
	while ((ch = getc(fp)) != EOF)
		{
		if ((ch == '\n') || (ch == '\r'))
			{
			if (ch == '\n')
				term = LINETERM_UNIX;
			else /* ch = '\r' */
				{
				/* First test below handles one-line MAC file */
				if (((nextCh = getc(fp)) == EOF) || (nextCh != '\n'))
					term = LINETERM_MAC;
				else
					term = LINETERM_DOS;
				}
			break;
			}
		}
	(void)fseek(fp, 0L, 0);		/* rewind */
	
	return (term);

}





int LongestLine (FILE *fp)

{

	int			ch, lineLength, longest;
	
	longest = 0;
	lineLength = 0;
	while ((ch = fgetc(fp)) != EOF)
		{
		if ((ch == '\n') || (ch == '\r'))
			{
			if (lineLength > longest)
				longest = lineLength;
			lineLength = 0;
			}
		else
			lineLength++;
		}
	rewind (fp);		/* rewind */
	
	return (longest);

}





void MrBayesPrint (char *format, ...)

{
	va_list ptr;

#	if defined (MPI_ENABLED)
	if (proc_id == 0)
		{
		if (echoMB == YES)
			{
			va_start (ptr, format);
			vprintf (format, ptr);
			va_end(ptr);
			fflush (stdout);
			}
		if (logToFile == YES)
			{
			if (logFileFp == NULL)
				printf ("%s   Could not print log output to file\n", spacer);
			else
				{
				va_start (ptr, format);
				vfprintf (logFileFp, format, ptr);
				va_end(ptr);
				fflush (logFileFp);
				}
			}
		}
#	else
	if (chainParams.redirect == NO)
		{
		if (echoMB == YES)
			{
			va_start (ptr, format);
			vprintf (format, ptr);
			va_end(ptr);
			fflush (stdout);
			}
		if (logToFile == YES)
			{
			if (logFileFp == NULL)
				{
				printf ("%s   Could not print log output to file\n", spacer);
				logToFile = NO;
				}
			else
				{
				va_start (ptr, format);
				vfprintf (logFileFp, format, ptr);
				va_end(ptr);
				fflush (logFileFp);
				}
			}
		}
#	endif
}




void MrBayesPrintf (FILE *f, char *format, ...)

{
	va_list                 ptr;

#	if defined (MPI_ENABLED)
	if (proc_id == 0)
		{
		va_start (ptr, format);
		vfprintf (f, format, ptr);
		va_end(ptr);
		fflush(f);
		}
#	else
	va_start (ptr, format);
	vfprintf (f, format, ptr);
	va_end(ptr);
	fflush(f);
#	endif
}





FILE *OpenBinaryFileR (char *name)

{

	FILE		*fp;

	if ((fp = fopen (name, "rb")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, name);
		return (NULL);
		}
	else
		return (fp);
	
}





FILE *OpenTextFileR (char *name)

{

	FILE		*fp;

	if ((fp = fopen (name, "r")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, name);
		return (NULL);
		}
	else
		return (fp);
	
}





FILE *OpenTextFileA (char *name)

{

	FILE		*fp;

	if ((fp = fopen (name, "a+")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, name);
		return (NULL);
		}
	else
		return (fp);
	
}





FILE *OpenTextFileW (char *name)

{

	FILE		*fp;

	if ((fp = fopen (name, "w+")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, name);
		return (NULL);
		}
	else
		return (fp);
	
}




int SafeFclose(FILE **fp) {
	int retval=-1;
#if defined MPI_ENABLED
	if (proc_id == 0) {
#endif
	if( fp!=NULL && (*fp)!=NULL ) 
		retval=fclose(*fp);
	*fp = NULL;
#if defined MPI_ENABLED
	}
#endif
	return retval;	
}





void *SafeMalloc(size_t s) {
        void *ptr = malloc(s);
        if(ptr==NULL)
                return NULL;
        return memset(ptr,0,s);
}





/* SafeStrcat: Allocate or reallocate target to fit result; assumes ptr is NULL if not allocated */
char *SafeStrcat (char **target, const char *source)
{
    if (*target == NULL)
        *target = (char *) malloc ((strlen(source)+1)*sizeof(char));
    else
        *target = (char *) realloc ((void *)*target, (size_t)(strlen(source)+strlen(*target)+1)*sizeof(char));

    if (*target)
        strcat(*target,source);

    return (*target);
}





/* StripComments: Strip possibly nested comments from the string s.
	Example: s="[[comment]]"-> s="comment" */
void StripComments (char *s)
{
	char	*t;
	int		inComment;

	inComment = 0;
	for (t=s; *s != '\0'; s++)
		{
		if (inComment == 0)
			{
			if (*s == '[')
				inComment++;
			else
				*t++ = *s;
			}
		else
			{
			if (*s == ']')
				inComment--;
            else if (*s == '[')
                inComment++;
			}
 		}
    *t = '\0';
}





