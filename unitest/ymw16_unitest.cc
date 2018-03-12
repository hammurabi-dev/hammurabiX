///
/// unit-test for namespace toolkit
/// feel free to add more rational testing blocks
///
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cgs_units_file.h>
#include <ap_err.h>
#include <fereg.h>
#include <param.h>
using namespace std;

template<typename T>
bool compare(const T &,const T &,const double &);

template<typename T>
bool compare(const T &a,const T &b,const double &precision){
    return (fabs(a-b)<precision);
}

int main(void){
    cout<<endl<<"YMW16 regular free-electron field unitest"<<endl;
    
    // load template params
    unique_ptr<Param> par = unique_ptr<Param> (new Param("reference/ymw16_unitest.xml"));
    
    unique_ptr<Grid_fereg> test_grid = unique_ptr<Grid_fereg> (new Grid_fereg("reference/ymw16_unitest.xml"));
    unique_ptr<FEreg> test_field = unique_ptr<FEreg> (new FEreg_ymw16());
    test_field->write_grid(par.get(),test_grid.get());
    //grid->export_grid();
    
    // read reference file
    unique_ptr<Grid_fereg> ref_grid = unique_ptr<Grid_fereg> (new Grid_fereg("reference/ymw16_unitest.xml"));
    unique_ptr<FEreg> ref_field = unique_ptr<FEreg> (new FEreg_ymw16());
    ref_grid->import_grid();
    
    if(ref_grid->full_size!=test_grid->full_size){
        ap_err("YMW16 ...... fail");
        exit(1);
    }
    for(decltype(ref_grid->full_size) i=0;i!=ref_grid->full_size;++i){
        if(!compare(ref_grid->fe[i],test_grid->fe[i],1e-5)){
            cout<<"ref "<<ref_grid->fe[i]<<endl
            <<"test "<<test_grid->fe[i]<<endl;
            ap_err("YMW16 ...... fail");
            exit(1);
        }
    }
    
    // if all testing blocks pass
    cout<<"YMW16 ...... pass"<<endl<<endl;
    return EXIT_SUCCESS;
}
