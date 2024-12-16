void latticeMeasure(double V, int Nconfig);
void endMeasure(double T, double B0, double V0, double phi, int Nmeasurment);


double calculateLayerEnergy(int k, double V);
double calculateDiagonalEnergy(int k, double V);
double calculateTrotterEnergy(int k);

double calculateLayerZeemanEnergy(int k);
double calculateMagEnergy(int k, int ix, int iy, int itri, int ispin);

double calculateLayerMagX(int k);
double calculateMagX(int k, int ix, int iy, int itri, int ispin);

double calculateLayerMagY(int k);
double calculateMagY(int k, int ix, int iy, int itri, int ispin);

double calculateLayerMag(int k);
double calculateMag(int k, int ix, int iy, int itri, int ispin);

int topologicalSector0(int k);
int topologicalSector1(int k);
int topologicalSector2(int k);

void magnetizationHistogramMeasure();
void orderParameterHistogramMeasure();

void save_lattice();
void save_ref();

bool weAreInStarPhaseTopologicalSector(int k);
