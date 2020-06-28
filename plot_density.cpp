void plot_density(){
    
    FILE *myfile = fopen("density_0000", "r");
    
    int Nx = 32, Ny = 32, Nz = 32;
    double Lx = 1.0, Ly = 1.0, Lz = 1.0;
    int i, j, k;
    double rho;
    
    TH3F *hist = new TH3F("glvoxel", "glvoxel", Nx, 0, Lx, Ny, 0, Ly, Nz, 0, Lz);
    
    while(!feof(myfile)){
        fscanf(myfile, "%d\t%d\t%d\t%le", &i, &j, &k, &rho);
        hist->SetBinContent(i, j, k, rho);
    }
    
    gStyle->SetCanvasPreferGL(1);
    gStyle->SetPalette(1);
    
    hist->Draw("box2z");
    
    fclose(myfile);
}


