#ifndef GLOBALS_H
#define GLOBALS_H

#define _INPUT_COMMENT_ '#'
#define _INPUT_DECLARATION_ '='

#define _FASTA_COMMENT_ '>'

#define _INPUT_SEPARATOR_ '|'
#define _OUTPUT_SEPARATOR_ '|'

#define _ERROR_INPUT_UNREADABLE_ "Input file impossible to open"
#define _ERROR_FASTA_UNREADABLE_ "Fasta file impossible to open"
#define _ERROR_MUTATION_TARGET_UNFINDABLE_ "Did not find mutation target"
#define _ERROR_OUTPUT_BUFFER_ "Error: trying to add data of a step that has already been written to the result file."

#define _MIN_OUTPUT_PRECISION_ 2

#define _INPUT_KEY_GENERATIONS_ "GEN"
#define _INPUT_KEY_REPLICAS_ "REP"
#define _INPUT_KEY_MARKER_SITES_ "SITES"
#define _INPUT_KEY_MODE_ "MODE"
#define _INPUT_KEY_MIGRATION_MODEL_ "MIG_MODEL"
#define _INPUT_KEY_MIGRATION_MODE_ "MIG_MODE"
#define _INPUT_KEY_MIGRATION_RATES_ "MIG_RATES"
#define _INPUT_KEY_MUTATION_RATES_ "MUT"
#define _INPUT_KEY_MUTATION_KIMURA_ "MUT_KIMURA"
#define _INPUT_KEY_MUTATION_FELSENSTEIN_ "MUT_FELSENSTEIN"
#define _INPUT_KEY_SELECTION_RATES_ "SEL"

#define _PARAM_NONE_ 0
#define _PARAM_MUTATIONS_ 1
#define _PARAM_MIGRATION_ 2
#define _PARAM_SELECTION_ 4
#define _PARAM_BOTTLENECK_ 8

#define _DEFAULT_MUTATION_RATE_ 1E-6
#define _MUTATION_MODEL_NONE_ 0
#define _MIGRATION_MODEL_NONE_ 0

#define _MUTATION_MODEL_CANTOR_ 1
#define _MUTATION_MODEL_KIMURA_ 2
#define _MUTATION_MODEL_FELSENSTEIN_ 3

#define _BOTTLENECK_START_TIME_ 20
#define _BOTTLENECK_END_TIME_ 40
#define _BOTTLENECK_POPULATION_REDUCTION_ 2.0

#define _COMPLETE_GRAPH_ 0
#define _STAR_ 1
#define _RING_ 2

#define _INPUT_USER_ 0
#define _RANDOM_ 1
#define _DEFAULT_EXCESS_ 1E6
#define _MIGRATION_OUTPUT_SEPARATOR_ "  "

#define _MIGRATION_DETAILED_OUTPUT_ true

#define strToInt [](const std::string& s){ return std::stoi(s); }
#define strToDouble [](const std::string& s){ return std::stod(s); }

#endif
