#include<iostream>
#include<cstdlib>
#include<cmath>
#include<TStyle.h>
#include<TCanvas.h>
#include<TH3F.h>
#include<TTimer.h>
#include<TText.h>

using namespace std;

double dt = 0.1;
int N = 128;
int Nx = N, Ny = N, Nz = N;
double Lx = 1.0, Ly = 1.0, Lz = 1.0;
int i, j, k;
double rho;
char name[64], output_name[64], time_text[64];
int frame = 0;
double t = 0;
double phi = 0;
FILE *myfile;

TCanvas *c1 = new TCanvas("c1", "Density", 800, 600);
TH3F *hist = new TH3F("Density", "Density", Nx, 0, Lx, Ny, 0, Ly, Nz, 0, Lz);
TText *text = new TText(0.5, 0.5, "");

void plot_density(){
    
    gStyle->SetCanvasPreferGL(1);
    gStyle->SetPalette(1);
    
    TTimer *timer = new TTimer(0);
    timer->SetCommand("anim()");
    timer->TurnOn();
    
}

void anim(){
    hist->Reset();
    
    sprintf(name, "./output/density_%04d", frame);
    if((myfile = fopen(name, "r")) == NULL){
        frame = 0;
        sprintf(name, "./output/density_%04d", frame);
        myfile = fopen(name, "r");
        t = 0;
    }
    
    while(!feof(myfile)){
        fscanf(myfile, "%d\t%d\t%d\t%le", &i, &j, &k, &rho);
        hist->SetBinContent(i, j, k, rho);
    }
    
    hist->SetMaximum(90000);
    hist->SetMinimum(0);
    hist->Draw("box2z");
    sprintf(time_text, "Density: frame %d", frame);
    hist->SetTitle(time_text);
    //gPad->SetPhi(phi);
    gPad->Modified();
    gPad->Update();
    sprintf(output_name, "./output/density_%04d.gif", frame);
    c1->SaveAs(output_name);
    
    fclose(myfile);
    phi+=0.2;
    frame+=1;
    
}
