#ifndef KVEC_H
#define KVEC_H

#include <stdlib.h>

#define kvec_t(type) struct { size_t n, m; type *a; }

#endif
