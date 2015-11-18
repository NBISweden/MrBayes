#ifndef __SUMPT_H__
#define __SUMPT_H__

/* struct to hold info about a .p file */
typedef struct
    {
    int     longestLineLength;
    int     headerLine;
    int     firstParamLine;
    int     numRows;
    int     numColumns;
    } SumpFileInfo;

/* struct to hold info about a model probability */
typedef struct
    {
    int     index;
    double  prob;
    } ModelProb;

/* struct to hold a parameter sample, possibly from multiple files */
typedef struct
    {
    MrBFlt **values;
    } ParameterSample;

/* function declarations */
int     AllocateParameterSamples (ParameterSample **parameterSamples, int numRuns, int numRows, int numColumns);
int     DoSump (void);
int     DoSumpParm (char *parmName, char *tkn);
int     DoSumSs (void);
int     DoSumSsParm (char *parmName, char *tkn);
int     ExamineSumpFile (char *fileName, SumpFileInfo *fileInfo, char ***headerNames, int *nHeaders);
int     FindHeader (char *token, char **headerNames, int nHeaders, int *index);
void    FreeParameterSamples (ParameterSample *parameterSamples);
int     GetHeaders (char ***headerNames, char *headerLine, int *nHeaders);
int     PrintPlot (MrBFlt *xVals, MrBFlt *yVals, int nSamples);
int     ReadParamSamples (char *fileName, SumpFileInfo *fileInfo, ParameterSample *parameterSamples, int runNo);

int     DoCompareTree (void);
int     DoCompareTreeParm (char *parmName, char *tkn);
int     DoCompRefTree (void);
int     DoSumt (void);
int     DoSumtParm (char *parmName, char *tkn);
int     DoSumtTree (void);
int     DoSumtTreeParm (char *parmName, char *tkn);
void    ResetTranslateTable (void);
int     ShowConTree (FILE *fp, PolyTree *t, int screenWidth, int showSupport);
void    ShowParts (FILE *fp, BitsLong *p, int nTaxaToShow);

#endif  /* __SUMPT_H__ */
