#include "AFocalSurface.h"
#include "AGeoBezierPgon.h"
#include "AGeoWinstonCone2D.h"
#include "AGeoWinstonConePoly.h"
#include "ALens.h"
#include "AMirror.h"
#include "AOpticsManager.h"
#include "ARayShooter.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoMatrix.h"
#include "TLegend.h"


// define useful units
static const Double_t cm = AOpticsManager::cm();
static const Double_t mm = AOpticsManager::mm();
static const Double_t um = AOpticsManager::um();
static const Double_t nm = AOpticsManager::nm();
static const Double_t m = AOpticsManager::m();

TGraph* ConeTrace(Int_t mode, bool del) {
  // mode == 0: Winston cone (built with AGeoWinstonConePoly)
  // mode == 1: Bezier cone (built with AGeoBezierPgon)

  AOpticsManager* manager = new AOpticsManager("manager", "SC");

  // Make the world
  TGeoBBox* worldbox = new TGeoBBox("worldbox", 10 * cm, 10 * cm, 10 * cm);
  AOpticalComponent* world = new AOpticalComponent("world", worldbox);
  manager->SetTopVolume(world);

  const Double_t kRin = 26 * mm;
  const Double_t kRout = 12.65 * mm;
  const Double_t kHeight = 11.5 * mm;
  TGeoRotation* rot30 = new TGeoRotation("rot30", 30, 0, 0);
  rot30->RegisterYourself();
  TGeoRotation* rot60 = new TGeoRotation("rot60", 60, 0, 0);
  rot60->RegisterYourself();
  TGeoRotation* rot120 = new TGeoRotation("rot120", 120, 0, 0);
  rot120->RegisterYourself();

  AGeoWinstonCone2D* coneV =
      new AGeoWinstonCone2D("coneV", kRin, kRout, kRin * 1.01);
  TGeoPgon* pgon = new TGeoPgon("pgon", 0, 360, 50, 4);
  pgon->DefineSection(0, -kHeight * 0.9999+(kRout/2), 0, kRout * 1.1);
  pgon->DefineSection(1, -kHeight * 0.5+(kRout/2), 0, kRin * 0.9);
  pgon->DefineSection(2, -kHeight * 0.+(kRout/2), 0, kRin * 1.001);
  pgon->DefineSection(3, kHeight * 0.9999+(kRout/2), 0, kRin * 1.001);

  AGeoWinstonConePoly* hexV = new AGeoWinstonConePoly("hexV", kRin, kRout, 50);
  TGeoCompositeShape* coneComp1 = 0;

  if (mode == 0) {
    coneComp1 = new TGeoCompositeShape("coneComp1", "pgon:rot30 - hexV");
  } else if (mode == 1) { 
    AGeoBezierPgon* hexB =
        new AGeoBezierPgon("hexB", 0, 360, 50, 100, kRin, kRout, kHeight);
    hexB->SetControlPoints(0.60, 0.25, 0.87, 0.40);
    AGeoBezierPgon* hexBout = 
        new AGeoBezierPgon("hexBout", 0, 360, 50, 100, kRin+0.2, kRout+0.2, kHeight);
    hexBout->SetControlPoints(0.60, 0.25, 0.87, 0.40);
    coneComp1 = new TGeoCompositeShape("coneComp1", "hexBout:rot30 - hexB:rot30");
  } //  if
  AMirror* coneMirror = new AMirror("coneMirror", coneComp1);
  world->AddNode(coneMirror, 1);

//  TGeoPgon* pgonPMT = new TGeoPgon("pgonPMT", 0, 360, 50, 2);
//  pgonPMT->DefineSection(0, -kHeight - 0.01 * mm, 0, kRout * 1.01);
//  pgonPMT->DefineSection(1, -kHeight, 0, kRout * 1.01);
//  AFocalSurface* pmt = new AFocalSurface("pmt", pgonPMT);
//  world->AddNode(pmt, 1, rot30);




//  TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
//  TGeoMedium *med = new TGeoMedium("MED",1,mat);
  TGeoSphere *roundPMT = new TGeoSphere("roundPMT", 0,kRout,0,90,0,360);
  AFocalSurface* pmt = new AFocalSurface("pmt", roundPMT);
  world->AddNode(pmt, 1);


  manager->CloseGeometry();

  if (mode == 1) {
    TCanvas* can1 = new TCanvas("can1", "can1");
    manager->GetTopVolume()->Draw("ogl");
  } //  if

  TGraph* graAeff = new TGraph;

  for (Double_t deg = 0.; deg < 80.; deg += 0.1) {
    Double_t Aeff = 0.;
    for (Double_t phi = 0.; phi < 30.; phi += 0.3) {
      TGeoTranslation* raytr = new TGeoTranslation(
          "raytr", 50 * mm * TMath::Sin(deg * TMath::DegToRad()), 0,
          50 * mm * TMath::Cos(deg * TMath::DegToRad()));
      TGeoRotation* rayrot = new TGeoRotation("rayrot", 90 - phi, 180 + deg, 0);

      TVector3 dir(0, 0, 1);
      Double_t lambda =
          400 * nm;  // does not affect the results because we have no lens
      // 1 photon per 1 mm^2
      const int kN = 1000;
      const double kD = 100 * mm;
      ARayArray* array =
          ARayShooter::RandomSquare(lambda, kD, kN, rayrot, raytr, &dir);
      Double_t dA = kD * kD / kN;

      manager->TraceNonSequential(*array);
      TObjArray* focused = array->GetFocused();

      for (Int_t j = 0; j <= focused->GetLast(); j++) {
        ARay* ray = (ARay*)(*focused)[j];
        if (!ray) continue;

        // Calculate the effective area from the number of focused photons
        Aeff += dA;

        if (mode == 1 and TMath::Abs(deg - 20.) < 0.001 and j < 50 and
            phi == 0) {
          TPolyLine3D* pol = ray->MakePolyLine3D();
          pol->SetLineColor(3);
          pol->SetLineWidth(1);
          pol->Draw();
        }  // if
      }    // j

      delete array;
    }  // phi

    graAeff->SetPoint(graAeff->GetN(), deg,
                      Aeff / (2 * TMath::Sqrt(3) * kRin * kRin) /
                          TMath::Cos(deg * TMath::DegToRad()));
  }  // deg

  if (del) {
    delete manager;
    manager = 0;
  }  // if

  return graAeff;
}

void winstonbeziertest() {
  TGraph* win = ConeTrace(0, true); //   Winston cone
  TGraph* bez = ConeTrace(1, false);  // Bezier cone

  TGraph* gra[2] = {win, bez};
  TCanvas* can2 = new TCanvas("can2", "can2");


  for (int i = 0; i < 2; i++) {
    gra[i]->GetXaxis()->SetTitle("Incident Angle (deg)");
    gra[i]->GetYaxis()->SetTitle("Collection Efficiency (%)");
    gra[i]->SetMarkerSize(0.3);
    gra[0]->SetMarkerColor(kRed);
    gra[1]->SetMarkerColor(kBlue);
    gra[i]->SetMarkerStyle(20);
    gra[i]->Draw(i == 0 ? "ap" : "p same");
    gra[i]->GetXaxis()->SetRangeUser(0, 89.);
    gra[i]->GetYaxis()->SetRangeUser(0, 100.);
  }  // i

  TLegend* leg = new TLegend(0.15, 0.15, 0.6, 0.25);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(132);
  leg->AddEntry(win, "Winston Cone", "p");
  leg->AddEntry(bez, "Bezier Cone", "p");
  leg->Draw();
}
