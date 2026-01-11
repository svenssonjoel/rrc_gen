/* MIT License */

/* Copyright (c) 2026 Joel Svensson */

/* Permission is hereby granted, free of charge, to any person obtaining a copy */
/* of this software and associated documentation files (the "Software"), to deal */
/* in the Software without restriction, including without limitation the rights */
/* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
/* copies of the Software, and to permit persons to whom the Software is */
/* furnished to do so, subject to the following conditions: */

/* The above copyright notice and this permission notice shall be included in all */
/* copies or substantial portions of the Software. */

/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
/* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, */
/* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE */
/* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
/* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, */
/* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE */
/* SOFTWARE. */
                  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <stdbool.h>
#include <string.h>

#define GEN_HEADER    0x401
#define GEN_NORMALIZE 0x402

struct option options[] = {
  {"help", no_argument, NULL, 'h'},
  {"tsym", required_argument, NULL, 't'},
  {"ntaps", required_argument, NULL, 'n'},
  {"window", required_argument, NULL, 'w'}, // number of symbols to cover with the filter.
  {"alpha", required_argument, NULL, 'a'},
  {"header", no_argument, NULL, GEN_HEADER},
  {"name", required_argument, NULL, 'o'},
  {"normalize", no_argument, NULL, GEN_NORMALIZE},
  {0,0,0,0}
};

// Defaults 
static int gen_filter_symbols_per_second  = 1000; // default 1000sym/s
static int gen_filter_taps_per_symbol     = 8;
static int gen_filter_window              = 11;
static float gen_filter_alpha             = 0.3f;
static bool gen_header                    = false;
static bool gen_normalize                 = false;

// File generation
static char *gen_name                     = "out";
static char *csv_ending                   = ".csv";
static char *csv_taps_ending              = "_taps.csv";
static char *h_ending                     = ".h";

void parse_opts(int argc, char **argv) {
  int c;
  opterr = 1;
  int opt_index = 0;
  while ((c = getopt_long(argc, argv, "hb:n:a:w:o:t:",options, &opt_index)) != -1) {
    switch(c) {
    case 'h':
      printf("Usage: %s [OPTION...]\n\n", argv[0]);
      printf("    -h           Prints this help\n");
      printf("    -t value     Set the symbols per seconds to value\n");
      printf("    -n value     Set number of taps per symbol to value\n");
      printf("    -a value     Set alpha value to value\n");
      printf("    -o filename  Write output to file filename\n");
      printf("\n");
      printf("    --header     Generate a C header containing the filter def\n");
      printf("    --normalize  Normalize filter so that total filter energy = 1\n");
             
      break;
    case 't':
      gen_filter_symbols_per_second = (int)atoi((char *)optarg);
      break;
    case 'a':
      gen_filter_alpha = strtof(optarg, NULL);
      break;
    case 'n':
      gen_filter_taps_per_symbol = (int)atoi((char *)optarg);
      break;
    case 'w':
      gen_filter_window = (int)atoi((char *)optarg);
      break;
    case 'o':
      gen_name = (char *)optarg;
      break;
    case GEN_HEADER:
      gen_header = true;
      break;
    case GEN_NORMALIZE:
      gen_normalize = true;
      break;
    }
  }
}


// calculate the filter's total "energy"
float total_energy(float *filter, int filter_len) {
  float energy = 0.0f;
  for (int i = 0; i < filter_len; i ++) { 
    energy += filter[i] * filter[i];
  }
  return energy;
}

// Normalization.
//
//   Keep in mind that mag = sqrt(energy) for energy normalization.
void normalize_filter(float *filter, int filter_len, float mag) {
  for (int i = 0; i < filter_len; i ++) {
    filter[i] = filter[i] / mag;
  }
}

/**  compute a root raised Cosine filter tap at time t 
 *   with Symbol time t_sym.  Alpha controls the rolloff.
 *  
 * \param t_sym Symbol time.
 * \param t current time.
 * \param alpha rolloff characteristic. 
 * \return filter tap value at time t. 
 */
float rrc(float t_sym, float t, float alpha) {
  float T = 1.0f / t_sym;
  float epsilon = 1e-10f;  // Small value to handle numerical precision

  // Special case: t = 0
  if (fabsf(t) < epsilon) {
    return (1.0f / sqrtf(T)) * (1.0f + alpha * (4.0f / M_PI - 1.0f));
  }

  // Special case: t = ±T/(4*alpha)
  if (fabsf(fabsf(t) - T / (4.0f * alpha)) < epsilon) {
    return (alpha / sqrtf(2.0f * T)) *
           ((1.0f + 2.0f / M_PI) * sinf(M_PI / (4.0f * alpha)) +
            (1.0f - 2.0f / M_PI) * cosf(M_PI / (4.0f * alpha)));
  }

  // General case
  float numerator = sinf(M_PI * t / T * (1.0f - alpha)) +
                    4.0f * alpha * t / T * cosf(M_PI * t / T * (1.0f + alpha));
  float denominator = M_PI * t / T * (1.0f - powf(4.0f * alpha * t / T, 2.0f));

  return (1.0f / sqrtf(T)) * (numerator / denominator);
}

int main(int argc, char **argv) {


  parse_opts(argc,argv);

  // Filter parameters
  float t_sym = gen_filter_symbols_per_second;
  float alpha = gen_filter_alpha;
  float T = 1.0f / t_sym;  // Symbol period
  
  // time of the filter is tf = T * gen_filter_window.
  // The number of taps to generate is N = gen_filter_window * gen_filter_taps_per_symbol + 1

  // For visualization of this filter
  // sample from t = - (tf/2) to t = tf/2
  // and we are a bit free to select a number of samples to grab for visualization purpose.
  // that is dt is chosen freely. can we scale it nicely?
    
  // For header generation
  // Sample from t = - (tf/2) to t = tf/2
  // at dt = tf/(N-1)
  // collect N samples.
    
  // Sample over 6 symbol periods (-3T to +3T)

  float tf = T * gen_filter_window;
  float t_start = -(tf/2);
  float t_end   = tf/2;

  // N is used for tap generation.
  int N = gen_filter_window * gen_filter_taps_per_symbol + 1;

  if (gen_header) {

    char filename[256];
    memset(filename,0, 256);

    FILE *fp = NULL;
      
    snprintf(filename,256, "%s%s", gen_name, h_ending);
     
    printf("Generating Header file for use with C code\n\n");

    float dt = T / gen_filter_taps_per_symbol;

    float *coeffs = (float *)malloc(N * sizeof(float));
    // TODO: Error handling here... 
    
    for (int i = 0; i < N; i++) {
      float t = t_start + i * dt;
      coeffs[i] = rrc(t_sym, t, alpha);
    }

    float energy = total_energy(coeffs, N);
  
    // This "energy" is a dimensionless quantity... 
    if (gen_normalize) {
      printf("Normalizing filter\n");
      printf("Total energy: %.6f\n", energy);
      printf("Normalization factor: %.6f\n", sqrtf(energy));
      normalize_filter(coeffs, N, sqrtf(energy));
      printf("Total energy normalized: %.6f\n", total_energy(coeffs, N));
    }

    fp = fopen(filename, "w");
    if (fp) {
      fprintf(fp, "// RRC Filter Coefficients (Energy Normalized)\n");
      fprintf(fp, "// Symbol rate: %d symbols/sec\n", gen_filter_symbols_per_second);
      fprintf(fp, "// Taps per symbol: %d\n", gen_filter_taps_per_symbol);
      fprintf(fp, "// Window: %d symbols\n", gen_filter_window);
      fprintf(fp, "// Alpha: %.3f\n", alpha);
      fprintf(fp, "// Number of taps: %d\n", N);
      fprintf(fp, "// Energy: %.6f\n\n", total_energy(coeffs, N));
      
      fprintf(fp, "int %s_len = %d;\n", gen_name, N);
      fprintf(fp, "float %s[%d] = {\n", gen_name, N);
      
      for (int i = 0; i < N; i++) {
        fprintf(fp, "  %.10ef%s\n", coeffs[i], (i < N-1) ? "," : "");
      }
      
      fprintf(fp, "};\n");
      fclose(fp);
    } else {
      printf("Failed to open file %s for writing\n", filename);
    }
    free(coeffs);

    printf("RRC filter header written to %s\n", filename);
    printf("Parameters: %d taps, alpha=%.2f\n", N, alpha);
  } else {

    char csvname[256];
    memset(csvname,0, 256);
    char csvtapsname[256];
    memset(csvtapsname,0, 256);

    FILE *fp_csv = NULL;
    FILE *fp_csvtaps = NULL;
      
    snprintf(csvname,256, "%s%s", gen_name, csv_ending);
    snprintf(csvtapsname,256, "%s%s", gen_name, csv_taps_ending);

    printf("Generating CSV for visualisation\n\n");

    // Generate smooth curve for visualization
    int num_samples = 1000; // arbitrary to make it look ok in vis
    float dt_vis = (t_end - t_start) / (num_samples - 1);
    float dt_tap = T / gen_filter_taps_per_symbol;

    // TODO error handling here...
    float *vis_coeffs = (float *)malloc(num_samples * sizeof(float));
    float *tap_coeffs = (float *)malloc(N * sizeof(float));

    // Calculate coeffs for visualisation.
    for (int i = 0; i < num_samples; i++) {
      float t = t_start + i * dt_vis;
      float value = rrc(t_sym, t, alpha);
      vis_coeffs[i] = value;
    }

    // Calculate the filter taps
    for (int i = 0; i < N; i++) {
      float t = t_start + i * dt_tap;
      tap_coeffs[i] = rrc(t_sym, t, alpha);
    }

    // If normalize, normalize both taps and vis using the
    // normalization factor calculated from the taps.
    // This way the visualisation makes sense relative the taps.
    if (gen_normalize) {
      float energy = total_energy(tap_coeffs, N);
      printf("Normalizing filter visualization\n");
      printf("  Total energy: %.6f\n", energy);
      printf("  Normalization factor: %.6f\n", sqrtf(energy));
      normalize_filter(vis_coeffs, num_samples, sqrtf(energy));
      printf("  Total energy normalized: %.6f\n", total_energy(vis_coeffs, num_samples));

      // Vis total energy will be larger... 
      
      printf("\nNormalizing filter taps\n");
      printf("  Total energy: %.6f\n", energy);
      printf("  Normalization factor: %.6f\n", sqrtf(energy));
      normalize_filter(tap_coeffs, N, sqrtf(energy));
      printf("  Total energy normalized: %.6f\n", total_energy(tap_coeffs, N));
    }
    
    fp_csv = fopen(csvname, "w");
    
    if (fp_csv) {
    
      fprintf(fp_csv, "t,rrc_value\n");

      for (int i = 0; i < num_samples; i++) {
        float t = t_start + i * dt_vis;
        float value = vis_coeffs[i];
        fprintf(fp_csv, "%.10e,%.10e\n", t, value);
      }
    } else {
      printf("Failed to open file %s for writing\n", csvname);
    }

    fp_csvtaps = fopen(csvtapsname, "w");

    if (fp_csvtaps) {
      fprintf(fp_csvtaps, "t, tap_value\n");

      float dt = T / gen_filter_taps_per_symbol;

      for (int i = 0; i < N; i++) {
        float t = t_start + i * dt_tap;
        fprintf(fp_csvtaps, "%.10e,%.10e\n", t, tap_coeffs[i]);
      }
      fclose(fp_csvtaps);
    } else {
      printf("Failed to open file %s for writing\n", csvtapsname); 
    }

    printf("Parameters: t_sym=%.1f, alpha=%.2f\n", t_sym, alpha);
    printf("Time range: %.6f to %.6f seconds\n", t_start, t_end);

    free(tap_coeffs);
    free(vis_coeffs);
  } 
return 0;
}

