#ifndef IO_H
#define IO_H

void write_header(FILE* fp, int n, float ww);
void write_frame_data(FILE* fp, int n, float* x, int* c);

#endif /* IO_H */
