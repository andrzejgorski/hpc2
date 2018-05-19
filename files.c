#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "files.h"


int count_lines_in_file(char* filename) {
    int lines = 0;
    char ch;
    FILE* fp;
    fp = fopen(filename, "r");
    while(!feof(fp)) {
        ch = fgetc(fp);
        if(ch == '\n')
        {
            lines++;
        }
    }
    fclose(fp);
    return lines;
}
