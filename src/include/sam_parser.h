#ifndef SAMParser_H
#define SAMParser_H

#include <vector>
#include <string>
#include <string.h>
#include <stdlib.h>

#include "common.h"
#include "record.h"

class SAMParser {
	FILE *input;
    std::string fname;
    Record currentRecord;

public:
	SAMParser (const std::string &filename);
	~SAMParser (void);

public:
	std::string readComment (void);
	bool readNext ();
	bool hasNext (void);
	size_t fpos (void);

public:
	void parse (Record &line);
	std::string head (void);
	Record next (void) ;
};

#endif
