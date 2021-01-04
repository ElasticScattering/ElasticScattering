
#define CL_FAIL_CONDITION(err, msg)                                                                                                                                                          \
    if (err != CL_SUCCESS) {                                                                                                                                                                    \
        std::cout << "\x1B[31mCL Error [" << CLErrorString(err) << "]\x1B[0m: " << msg << "\n\tat: " << "L" << __LINE__ << ": " << __FILE__ << "\n\tin function: " << __FUNCTION__ << std::endl;   \
        exit(0);                                                                                                                                                                                \
    }

#define FAIL_CONDITION(cond, msg)                                                                                                                            \
	if (cond) {                                                                                                                                                 \
		std::cout << "\x1B[30mError\x1B[0m: " << msg << "\n\tat: " << "L" << __LINE__ << ": " << __FILE__ << "\n\tin function: " << __FUNCTION__ << std::endl;  \
		exit(0);                                                                                                                                                \
	}

#define CFG_EXIT_CONDITION(cond, msg)                                     \
	if (cond) {                                                           \
		std::cout << "\x1B[95mConfig Error\x1B[0m: " << msg << std::endl; \
		exit(0);                                                          \
	}