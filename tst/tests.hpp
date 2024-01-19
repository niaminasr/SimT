
#ifndef TST_TESTS
#define TST_TESTS

#include <iostream>
#include <cstdlib>
#include "geometry/Point.hpp"
#include "geometry/Form.hpp"

#include "geometry/Grid.hpp"
#include "solver/System.hpp"

unsigned long test_counter = 0;
unsigned long failed_counter = 0;
bool is_err = false;

#define TEST(boolean)   test_counter++;\
                        if(!(boolean)){\
                            failed_counter++;\
                            std::cout<<"[FAILED] test "<<test_counter<<" at line "<<__LINE__<<" !"<<std::endl;\
                        }else{\
                            std::cout<<"[PASSED] test "<<test_counter<<" at line "<<__LINE__<<std::endl;\
                        }

#define RETURN_TEST     std::cout<<"("<<test_counter-failed_counter<<"/"<<test_counter<<") tests PASSED"<<std::endl;\
                        return failed_counter;

void setErr(){
    is_err = true;
}

bool isErr(){
    bool result = is_err;
    is_err = false;
    return result;
}

#endif