# ROOTLibs

## Style
```c++
using namespace MGROOT;
MGROOT::LoadDefaultEnvironment();
Hist::AddDirectory();
```

## PdfEditor
```c++
using namespace MGROOT;
PdfEditor editor(Window(WindowSize::kSliceLR), "PdfPainter");
editor.create("Title", 2, 2);
editor.cd(1, PadAxis(0, 1));
editor.save();
```

## Hist
```c++
using namespace MGROOT;
Axis AXmom("Momentum [GeV]", 40, 0.5, 1000.0, AxisScale::kLog);
Hist* hist = new Hist("Hist", HistAxis(AXmom, "Events/Bin"));
hist->style();
hist->draw();


Hist::Load();
```
