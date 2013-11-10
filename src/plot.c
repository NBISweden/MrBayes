/*
 *  MrBayes 3
 *
 *  (c) 2002-2013
 *
 *  John P. Huelsenbeck
 *  Dept. Integrative Biology
 *  University of California, Berkeley
 *  Berkeley, CA 94720-3140
 *  johnh@berkeley.edu
 *
 *  Fredrik Ronquist
 *  Swedish Museum of Natural History
 *  Box 50007
 *  SE-10405 Stockholm, SWEDEN
 *  fredrik.ronquist@nrm.se
 *
 *  With important contributions by
 *
 *  Paul van der Mark (paulvdm@sc.fsu.edu)
 *  Maxim Teslenko (maxim.teslenko@nrm.se)
 *
 *  and by many users (run 'acknowledgements' to see more info)
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

#include "globals.h"
#include "command.h"
#include "plot.h"
#include "sump.h"
#include "utils.h"
#if defined(__MWERKS__)
#include "SIOUX.h"
#endif

const char* const svnRevisionPlotC="$Rev$";   /* Revision keyword which is expended/updated by svn on each commit/update*/

/* local (to this file) */
int     foundCurly;
char    *plotTokenP;/* plotToken[CMD_STRING_LENGTH];*/






int DoPlot (void)

{


    int             i, n, nHeaders, burnin, len, longestHeader, whichIsX, whichIsY, numPlotted;
    char            temp[100], **headerNames = NULL;
    SumpFileInfo    fileInfo;
    ParameterSample *parameterSamples;
    
#   if defined (MPI_ENABLED)
    if (proc_id != 0)
        return NO_ERROR;
#   endif

    /* initialize values */
    headerNames = NULL;
    nHeaders = 0;
    parameterSamples = NULL;

    /* tell user we are ready to go */
    MrBayesPrint ("%s   Plotting parameters in file %s ...\n", spacer, plotParams.plotFileName);
    
    /* examine plot file */
    if (ExamineSumpFile (plotParams.plotFileName, &fileInfo, &headerNames, &nHeaders) == ERROR)
        return ERROR;
        
    /* Calculate burn in */
    burnin = fileInfo.firstParamLine - fileInfo.headerLine - 1;
        
    /* tell the user that everything is fine */
    MrBayesPrint ("%s   Found %d parameter lines in file \"%s\"\n", spacer, fileInfo.numRows + burnin, plotParams.plotFileName);
    if (burnin > 0)
        MrBayesPrint ("%s   Of the %d lines, %d of them will be summarized (starting at line %d)\n", spacer, fileInfo.numRows+burnin, fileInfo.numRows, fileInfo.firstParamLine);
    else
        MrBayesPrint ("%s   All %d lines will be summarized (starting at line %d)\n", spacer, fileInfo.numRows, fileInfo.firstParamLine);
    MrBayesPrint ("%s   (Only the last set of lines will be read, in case multiple\n", spacer);
    MrBayesPrint ("%s   parameter blocks are present in the same file.)\n", spacer);
        
    /* allocate space to hold parameter information */
    if (AllocateParameterSamples (&parameterSamples, 1, fileInfo.numRows, fileInfo.numColumns) == ERROR)
        goto errorExit;

    /* Now we read the file for real. First, rewind file pointer to beginning of file... */
    if (ReadParamSamples (plotParams.plotFileName, &fileInfo, parameterSamples, 0) == ERROR)
        goto errorExit;
                    
    /* get length of longest header */
    longestHeader = 9; /* 9 is the length of the word "parameter" (for printing table) */
    for (i=0; i<nHeaders; i++)
        {
        len = (int) strlen(headerNames[i]);
        if (len > longestHeader)
            longestHeader = len;
        }
        
    /* print x-y plot of parameter vs. generation */
    whichIsX = -1;
    for (i=0; i<nHeaders; i++)
        {
        if (IsSame (headerNames[i], "Gen") == SAME)
            whichIsX = i;
        }       
        
    if (whichIsX < 0)
        {
        MrBayesPrint ("%s   Could not find a column labelled \"Gen\" \n", spacer);
        goto errorExit;
        }
        
    numPlotted = 0;
    for (n=0; n<nHeaders; n++)
        {
        strcpy (temp, headerNames[n]);
        whichIsY = -1;
        if (!strcmp(plotParams.match, "Perfect"))
            {
            if (IsSame (temp, plotParams.parameter) == SAME)
                whichIsY = n;
            }
        else if (!strcmp(plotParams.match, "All"))
            {
            whichIsY = n;
            }
        else
            {
            if (IsSame (temp, plotParams.parameter) == CONSISTENT_WITH)
                whichIsY = n;
            }
            
        if (whichIsY >= 0 && whichIsX != whichIsY)
            {
            MrBayesPrint ("\n%s   Rough trace plot of parameter %s:\n", spacer, headerNames[whichIsY]);
            if (PrintPlot (parameterSamples[whichIsX].values[0], parameterSamples[whichIsY].values[0], fileInfo.numRows) == ERROR)
                goto errorExit;
            numPlotted++;
            }
        }
        
    if (numPlotted == 0)
        {
        MrBayesPrint ("%s   Did not find any parameters matching \"%s\" to plot\n", spacer, plotParams.parameter);
        }
            
    /* free memory */
    for (i=0; i<nHeaders; i++)
        free (headerNames[i]);
    free(headerNames);
    FreeParameterSamples(parameterSamples);

    expecting = Expecting(COMMAND);

    return (NO_ERROR);
    
errorExit:

    /* free memory */
    for (i=0; i<nHeaders; i++)
        free (headerNames[i]);
    free(headerNames);
    FreeParameterSamples(parameterSamples);

    expecting = Expecting(COMMAND);

    return (ERROR);
}





int DoPlotParm (char *parmName, char *tkn)

{

    int         tempI;
    MrBFlt      tempD;
    char        tempStr[100];

    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A matrix must be specified before sumt can be used\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            expecting  = Expecting(PARAMETER);
            expecting |= Expecting(SEMICOLON);
            }
        /* set Filename (plotParams.plotFileName) ***************************************************/
        else if (!strcmp(parmName, "Filename"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                strcpy (plotParams.plotFileName, tkn);
                MrBayesPrint ("%s   Setting plot filename to %s\n", spacer, plotParams.plotFileName);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Relburnin (chainParams.relativeBurnin) ********************************************************/
        else if (!strcmp(parmName, "Relburnin"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.relativeBurnin = YES;
                    else
                        chainParams.relativeBurnin = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Relburnin\n", spacer);
                    return (ERROR);
                    }
                if (chainParams.relativeBurnin == YES)
                    MrBayesPrint ("%s   Using relative burnin (a fraction of samples discarded).\n", spacer);
                else
                    MrBayesPrint ("%s   Using absolute burnin (a fixed number of samples discarded).\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                return (ERROR);
                }
            }
        /* set Burnin (chainParams.chainBurnIn) *******************************************************/
        else if (!strcmp(parmName, "Burnin"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                chainParams.chainBurnIn = tempI;
                MrBayesPrint ("%s   Setting burnin to %d\n", spacer, chainParams.chainBurnIn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Burninfrac (chainParams.burninFraction) ************************************************************/
        else if (!strcmp(parmName, "Burninfrac"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                if (tempD < 0.01)
                    {
                    MrBayesPrint ("%s   Burnin fraction too low (< 0.01)\n", spacer);
                    return (ERROR);
                    }
                if (tempD > 0.50)
                    {
                    MrBayesPrint ("%s   Burnin fraction too high (> 0.50)\n", spacer);
                    return (ERROR);
                    }
                chainParams.burninFraction = tempD;
                MrBayesPrint ("%s   Setting burnin fraction to %.2f\n", spacer, chainParams.burninFraction);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                {
                return (ERROR);
                }
            }
        /* set Parameter (plotParams.parameter) *******************************************************/
        else if (!strcmp(parmName, "Parameter"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                strcpy (plotParams.parameter, tkn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Parameter (plotParams.match) *******************************************************/
        else if (!strcmp(parmName, "Match"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    strcpy (plotParams.match, tempStr);
                else
                    return (ERROR);
                
                MrBayesPrint ("%s   Setting plot matching to %s\n", spacer, plotParams.match);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }

        else
            return (ERROR);
        }

    return (NO_ERROR);

}

