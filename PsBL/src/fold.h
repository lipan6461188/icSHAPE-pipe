#ifndef FOLD_H
#define FOLD_H

#include "pan_type.h"
#include "exceptions.h"
#include <HybridRNA.h>

namespace pan
{

sp<StringArray> get_dot_bracket(RNA &ct);
sp<StringArray> fold_two_seq(const string &seq_1, const string &seq_2, const bool forbidIntramolecular = false);

/* number starts from 1, returned pairs is 1-base coorinates */
vector< pair<uLONG, uLONG> > get_pairs(RNA &ct, const uLONG number = 1);


}
#endif