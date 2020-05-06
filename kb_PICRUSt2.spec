/*
A KBase module: kb_PICRUSt2
*/

module kb_PICRUSt2 {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_picrust2_pipeline(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
