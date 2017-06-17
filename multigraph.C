//Create and Draw a TMultiGraph
//Author:: Rene Brun
{
   gStyle->SetOptFit();
   TCanvas *c1 = new TCanvas("c1","multigraph",700,500);
   c1->SetGrid();

      // draw a frame to define the range
   TMultiGraph *mg = new TMultiGraph();

      // create first graph
   const Int_t n1 = 3;
   Double_t px1[] = {0.5 , 1.5 , 2.5};
   Double_t py1[] = {187.299 , 0. , 4.20218};
   Double_t ex1[] = {0.5 , 0.5 , 0.5};
   Double_t ey1[] = {34.4818 , 0. , 2.97367};
   TGraphErrors *gr1 = new TGraphErrors(n1,px1,py1,ex1,ey1);
   gr1->SetTitle("MC Prediction");
   gr1->SetFillStyle(0);
   gr1->SetMarkerColor(kBlue);
   gr1->SetLineColor(kBlue);
   gr1->SetMarkerStyle(21);
   gr1->SetLineWidth(2);
   mg->Add(gr1);


      // create second graph
   const Int_t n2 = 3;
   Float_t x2[]  = {0.5 , 1.5 , 2.5};
   Float_t y2[]  = {315.695 , -1.75529 , 42.6024};
   Float_t ex2[] = {0.5 , 0.5 , 0.5};
   Float_t ey2[] = {107.362 , 2.87756 , 13.8978};
   TGraphErrors *gr2 = new TGraphErrors(n2,x2,y2,ex2,ey2);
   gr2->SetTitle("Fake Estimation - Closure Test");
   gr2->SetFillStyle(0);
   gr2->SetMarkerColor(kRed);
   gr2->SetLineColor(kRed);
   gr2->SetMarkerStyle(20);
   gr2->SetLineWidth(3);
   gr2->SetLineStyle(3);
   mg->Add(gr2);

   mg->Draw("ap");

     //force drawing of canvas to generate the fit TPaveStats
   //c1->Update();
   c1->BuildLegend();
   c1->SetLogy();
   return c1;
}
