/* typedefs */

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
    }   ModelProb;

/* struct to hold a parameter sample, possibly from multiple files */
typedef struct {
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


