#ifndef UTILS_TYPES____
#define UTILS_TYPES____

#include <chrono>
#include <climits>
#include <functional>
#include <stdlib.h>

#define NOT_EXIST UINT_MAX
#define UNMATCHED UINT_MAX
#define NOFAILNODE -1

// #define FullCoverage

// #define CP2LE



#define ALL
#define DSQL1 
//#define DSQL2
#define DSQL3

typedef unsigned int uint;
typedef uint32_t VertexID;
typedef uint32_t LabelID;
#define INVALID_VERTEX_ID 100000000

#define LOG(message) std::cout << message << std::endl;


#endif //UTILS_TYPES____
