; Foci Analyzer: Windows setup (Inno Setup 6). Build with installer\build-windows.ps1, which collects the
; files in installer\build\payload and passes the version.
;
; The setup holds the program and uv. After copying them it runs installer\install.py (through
; run-install.cmd), which downloads Python and the packages, and the networks the user chose. Everything
; ends up in one folder per user ({localappdata}\FociAnalyzer, no administrator rights); uninstalling
; removes that folder. The user's settings (%USERPROFILE%\.foci_analyzer) stay.
;
; Unattended (for example by IT):  FociAnalyzer-Setup.exe /VERYSILENT /NETWORKS=stardist,spotiflow
;   /NETWORKS=all|none|a comma-separated list; /OWN="cellpose=C:\envs\cp\Scripts\python.exe" (one network).

#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif

[Setup]
AppId={{6C1F7E52-3A0B-4C55-9E4D-F0C1A11A1A11}
AppName=Foci Analyzer
AppVersion={#AppVersion}
AppPublisher=Netherlands Cancer Institute - BioImaging Facility
AppPublisherURL=https://github.com/BioImaging-NKI
DefaultDirName={localappdata}\FociAnalyzer
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
OutputDir=build
OutputBaseFilename=FociAnalyzer-Setup-{#AppVersion}
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
SetupIconFile=FociAnalyzer.ico
WizardImageFile=wizard-large-100.bmp,wizard-large-125.bmp,wizard-large-150.bmp,wizard-large-200.bmp
WizardSmallImageFile=wizard-small-100.bmp,wizard-small-125.bmp,wizard-small-150.bmp,wizard-small-200.bmp
UninstallDisplayIcon={app}\program\installer\FociAnalyzer.ico
UninstallDisplayName=Foci Analyzer
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible

[Tasks]
Name: "desktopicon"; Description: "Put a shortcut on the desktop"; GroupDescription: "Shortcuts:"

[Files]
Source: "build\payload\*"; DestDir: "{app}\program"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "build\uv.exe"; DestDir: "{app}\bin"; Flags: ignoreversion

[Icons]
Name: "{userprograms}\Foci Analyzer"; Filename: "{app}\Foci Analyzer.exe"; WorkingDir: "{app}"; IconFilename: "{app}\program\installer\FociAnalyzer.ico"; Comment: "Foci in nuclei and cells"
Name: "{userdesktop}\Foci Analyzer"; Filename: "{app}\Foci Analyzer.exe"; WorkingDir: "{app}"; IconFilename: "{app}\program\installer\FociAnalyzer.ico"; Tasks: desktopicon

[UninstallDelete]
; the environments, Python, the download cache and the log were made by install.py, not copied by the setup
Type: filesandordirs; Name: "{app}"

[Run]
Filename: "{app}\Foci Analyzer.exe"; Description: "Start Foci Analyzer"; Flags: postinstall nowait skipifsilent unchecked

[Code]
const
  NetCount = 6;
  { disk space in tenths of a GB once installed (the same numbers as SIZE in install.py) }
  SizeProgram = 5;
  SizeTorchGpu = 44;
  SizeTorchCpu = 10;
  EM_SCROLLCARET = $00B7;
var
  NetPage: TInputOptionWizardPage;
  TotalLabel: TNewStaticText;
  OwnPage: TWizardPage;
  OwnIntro: TNewStaticText;
  OwnLabel: array[0..5] of TNewStaticText;
  OwnEdit: array[0..5] of TNewEdit;
  OwnButton: array[0..5] of TNewButton;
  Names, Shown: array[0..5] of String;
  Sizes: array[0..5] of Integer;
  HasGpu: Boolean;
  LogMemo: TNewMemo;
  ProgressFile: String;
  LinesShown: Integer;

function SetTimer(hWnd, nIDEvent, uElapse, lpTimerFunc: Longword): Longword;
external 'SetTimer@user32.dll stdcall';
function KillTimer(hWnd, nIDEvent: Longword): Bool;
external 'KillTimer@user32.dll stdcall';

function GB(Tenths: Integer): String;
begin
  Result := IntToStr(Tenths div 10) + '.' + IntToStr(Tenths mod 10) + ' GB';
end;

function TorchSize: Integer;
begin
  if HasGpu then Result := SizeTorchGpu else Result := SizeTorchCpu;
end;

{ Cellpose, InstanSeg, Spotiflow and Trackastra (1 to 4) share one PyTorch; Cellpose 4 (5) has its own }
procedure NetsChanged(Sender: TObject);
var
  i, Total: Integer;
  Torch: Boolean;
begin
  Total := SizeProgram;
  Torch := False;
  for i := 0 to NetCount - 1 do
    if NetPage.Values[i] then
    begin
      Total := Total + Sizes[i];
      if (i >= 1) and (i <= 4) then Torch := True;
      if i = 5 then Total := Total + TorchSize;
    end;
  if Torch then Total := Total + TorchSize;
  TotalLabel.Caption := 'Together with the program: about ' + GB(Total) + ' on disk. Some models are downloaded ' +
    'when they are first used (CellposeSAM: 1.2 GB).';
end;

procedure BrowseClick(Sender: TObject);
var
  i: Integer;
  FileName: String;
begin
  i := TNewButton(Sender).Tag;
  FileName := OwnEdit[i].Text;
  if GetOpenFileName('The python.exe of your ' + Shown[i] + ' environment', FileName, '',
                     'python.exe|python.exe|All files|*.*', 'exe') then
    OwnEdit[i].Text := FileName;
end;

procedure InitializeWizard;
var
  i: Integer;
  Torch: String;
begin
  Names[0] := 'stardist'; Names[1] := 'cellpose'; Names[2] := 'instanseg'; Names[3] := 'spotiflow'; Names[4] := 'trackastra'; Names[5] := 'cellpose4';
  Shown[0] := 'StarDist'; Shown[1] := 'Cellpose'; Shown[2] := 'InstanSeg'; Shown[3] := 'Spotiflow'; Shown[4] := 'Trackastra'; Shown[5] := 'Cellpose 4';
  Sizes[0] := 16; Sizes[1] := 1; Sizes[2] := 3; Sizes[3] := 2; Sizes[4] := 2; Sizes[5] := 1;
  HasGpu := FileExists(ExpandConstant('{sys}\nvcuda.dll')) or FileExists(ExpandConstant('{sys}\nvidia-smi.exe'));
  if HasGpu then
    Torch := 'PyTorch for the NVIDIA graphics card (' + GB(TorchSize) + ')'
  else
    Torch := 'PyTorch for the processor (' + GB(TorchSize) + '; no NVIDIA graphics card found)';

  NetPage := CreateInputOptionPage(wpSelectTasks, 'Neural networks',
    'Which networks should be installed?',
    'Foci Analyzer works without them (threshold segmentation, watershed foci), but they find nuclei and foci ' +
    'better. Cellpose, InstanSeg, Spotiflow and Trackastra share one ' + Torch + '; Cellpose 4 has its own. ' +
    'StarDist runs on the processor. Networks you leave unticked can be taken from a Python environment of your ' +
    'own on the next page, or added later.', False, False);
  NetPage.Add('StarDist: nuclei, fast and robust (recommended) - ' + GB(Sizes[0]));
  NetPage.Add('Cellpose: nuclei or whole cells, also in 3D - ' + GB(Sizes[1]) + ' + PyTorch');
  NetPage.Add('InstanSeg: nuclei or whole cells, fast; also made for tissue - ' + GB(Sizes[2]) + ' + PyTorch');
  NetPage.Add('Spotiflow: foci detection by a network - ' + GB(Sizes[3]) + ' + PyTorch');
  NetPage.Add('Trackastra: tracking of cells in time-lapses - ' + GB(Sizes[4]) + ' + PyTorch');
  NetPage.Add('Cellpose 4: CellposeSAM, the newest model; slow without a graphics card - ' + GB(Sizes[5] + TorchSize));
  for i := 0 to NetCount - 2 do
    NetPage.Values[i] := True;
  TotalLabel := TNewStaticText.Create(NetPage);
  TotalLabel.Parent := NetPage.Surface;
  TotalLabel.AutoSize := False;
  TotalLabel.WordWrap := True;
  TotalLabel.Left := 0;
  TotalLabel.Width := NetPage.SurfaceWidth;
  TotalLabel.Height := ScaleY(34);
  TotalLabel.Top := NetPage.SurfaceHeight - TotalLabel.Height;
  TotalLabel.Anchors := [akLeft, akRight, akBottom];
  TotalLabel.Font.Style := [fsBold];
  NetPage.CheckListBox.Height := TotalLabel.Top - NetPage.CheckListBox.Top - ScaleY(6);
  NetPage.CheckListBox.OnClickCheck := @NetsChanged;
  NetsChanged(nil);

  OwnPage := CreateCustomPage(NetPage.ID, 'Your own environments (optional)',
    'Use a network from a Python environment of your own?');
  OwnIntro := TNewStaticText.Create(OwnPage);
  OwnIntro.Parent := OwnPage.Surface;
  OwnIntro.AutoSize := False;
  OwnIntro.WordWrap := True;
  OwnIntro.Width := OwnPage.SurfaceWidth;
  OwnIntro.Height := ScaleY(84);
  OwnIntro.Anchors := [akLeft, akTop, akRight];
  OwnIntro.Caption := 'These networks are not installed. If one of them is installed in a Python environment of ' +
    'your own, give the python.exe of that environment (for example C:\PythonProjects\Cellpose4\.venv\Scripts\' +
    'python.exe): Foci Analyzer then uses it. Leave a box empty otherwise. This can be changed later on the ' +
    'Settings page of Foci Analyzer.';
  for i := 0 to NetCount - 1 do
  begin
    OwnLabel[i] := TNewStaticText.Create(OwnPage);
    OwnLabel[i].Parent := OwnPage.Surface;
    OwnLabel[i].Caption := Shown[i] + ':';
    OwnEdit[i] := TNewEdit.Create(OwnPage);
    OwnEdit[i].Parent := OwnPage.Surface;
    OwnEdit[i].Left := ScaleX(80);
    OwnEdit[i].Width := OwnPage.SurfaceWidth - ScaleX(80) - ScaleX(88);
    OwnEdit[i].Anchors := [akLeft, akTop, akRight];
    OwnButton[i] := TNewButton.Create(OwnPage);
    OwnButton[i].Parent := OwnPage.Surface;
    OwnButton[i].Caption := 'Browse...';
    OwnButton[i].Width := ScaleX(80);
    OwnButton[i].Height := OwnEdit[i].Height + ScaleY(2);
    OwnButton[i].Left := OwnPage.SurfaceWidth - OwnButton[i].Width;
    OwnButton[i].Anchors := [akTop, akRight];
    OwnButton[i].Tag := i;
    OwnButton[i].OnClick := @BrowseClick;
  end;

  { the installation's messages, shown under the progress bar instead of in a console window }
  LogMemo := TNewMemo.Create(WizardForm);
  LogMemo.Parent := WizardForm.InstallingPage;
  LogMemo.Left := WizardForm.ProgressGauge.Left;
  LogMemo.Width := WizardForm.ProgressGauge.Width;
  LogMemo.Top := WizardForm.ProgressGauge.Top + WizardForm.ProgressGauge.Height + ScaleY(12);
  LogMemo.Height := WizardForm.InstallingPage.Height - LogMemo.Top - ScaleY(4);
  LogMemo.Anchors := [akLeft, akTop, akRight, akBottom];
  LogMemo.ReadOnly := True;
  LogMemo.ScrollBars := ssVertical;
  LogMemo.WordWrap := False;
  LogMemo.Font.Name := 'Consolas';
  LogMemo.Font.Size := 8;
  LogMemo.Visible := False;
end;

{ only the networks that are not installed are asked for }
procedure CurPageChanged(CurPageID: Integer);
var
  i, Y: Integer;
begin
  if CurPageID = OwnPage.ID then
  begin
    Y := OwnIntro.Top + OwnIntro.Height + ScaleY(10);
    for i := 0 to NetCount - 1 do
    begin
      OwnLabel[i].Visible := not NetPage.Values[i];
      OwnEdit[i].Visible := not NetPage.Values[i];
      OwnButton[i].Visible := not NetPage.Values[i];
      if not NetPage.Values[i] then
      begin
        OwnEdit[i].Top := Y;
        OwnButton[i].Top := Y - ScaleY(1);
        OwnLabel[i].Top := Y + ScaleY(3);
        Y := Y + OwnEdit[i].Height + ScaleY(10);
      end;
    end;
  end;
end;

function ShouldSkipPage(PageID: Integer): Boolean;
var
  i: Integer;
begin
  Result := False;
  if PageID = OwnPage.ID then
  begin
    Result := True;
    for i := 0 to NetCount - 1 do
      if not NetPage.Values[i] then Result := False;
  end;
end;

function OwnPython(i: Integer): String;
begin
  Result := '';
  if not NetPage.Values[i] then Result := Trim(OwnEdit[i].Text);
end;

function Quote(S: String): String;
begin
  Result := '"' + S + '"';
end;

function InstallArguments: String;
var
  i: Integer;
  Nets: String;
begin
  if ExpandConstant('{param:networks|}') <> '' then
  begin
    Result := '--networks ' + ExpandConstant('{param:networks|}');
    if ExpandConstant('{param:own|}') <> '' then
      Result := Result + ' --own ' + Quote(ExpandConstant('{param:own|}'));
    Exit;
  end;
  Nets := '';
  for i := 0 to NetCount - 1 do
    if NetPage.Values[i] then
    begin
      if Nets <> '' then Nets := Nets + ',';
      Nets := Nets + Names[i];
    end;
  if Nets = '' then Nets := 'none';
  Result := '--networks ' + Nets;
  for i := 0 to NetCount - 1 do
    if OwnPython(i) <> '' then
      Result := Result + ' --own ' + Quote(Names[i] + '=' + OwnPython(i));
end;

function NextButtonClick(CurPageID: Integer): Boolean;
var
  i: Integer;
begin
  Result := True;
  if CurPageID = OwnPage.ID then
    for i := 0 to NetCount - 1 do
      if (OwnPython(i) <> '') and not FileExists(OwnPython(i)) then
      begin
        MsgBox('There is no file ' + OwnPython(i) + ' (for ' + Shown[i] + ').', mbError, MB_OK);
        Result := False;
        Exit;
      end;
end;

{ install.py writes the progress file: per mille done, the step, the number of messages, the last messages }
procedure ShowProgress;
var
  Lines: TArrayOfString;
  i, Count, New, Last: Integer;
begin
  if not LoadStringsFromFile(ProgressFile, Lines) then Exit;
  Last := GetArrayLength(Lines) - 1;
  if Last < 2 then Exit;
  if StrToIntDef(Lines[0], -1) >= 0 then
    WizardForm.ProgressGauge.Position := 30 + StrToIntDef(Lines[0], 0) * 97 div 100;
  WizardForm.StatusLabel.Caption := Lines[1];
  Count := StrToIntDef(Lines[2], LinesShown);
  New := Count - LinesShown;
  if New > Last - 2 then New := Last - 2;
  if New > 0 then
  begin
    for i := Last - New + 1 to Last do
      LogMemo.Lines.Add(Lines[i]);
    while LogMemo.Lines.Count > 3000 do
      LogMemo.Lines.Delete(0);
    LogMemo.SelStart := Length(LogMemo.Text);
    SendMessage(LogMemo.Handle, EM_SCROLLCARET, 0, 0);
  end;
  LinesShown := Count;
end;

procedure ProgressTimer(Arg1, Arg2, Arg3, Arg4: Longword);
begin
  ShowProgress;
end;

procedure CurStepChanged(CurStep: TSetupStep);
var
  Code: Integer;
  Timer: Longword;
  Script: String;
begin
  if CurStep = ssPostInstall then
  begin
    ProgressFile := ExpandConstant('{app}\install-progress.txt');
    DeleteFile(ProgressFile);
    LinesShown := 0;
    WizardForm.ProgressGauge.Min := 0;
    WizardForm.ProgressGauge.Max := 1000;
    WizardForm.ProgressGauge.Position := 30;
    WizardForm.FilenameLabel.Caption := '';
    WizardForm.StatusLabel.Caption := 'Downloading Python, the program''s packages and the networks...';
    LogMemo.Visible := True;
    Script := ExpandConstant('{app}\program\installer\run-install.cmd');
    Timer := SetTimer(0, 0, 500, CreateCallback(@ProgressTimer));
    if not Exec(ExpandConstant('{cmd}'), '/c "' + Quote(Script) + ' ' + InstallArguments + ' --progress ' +
                Quote(ProgressFile) + '"', '', SW_HIDE, ewWaitUntilTerminated, Code) then
      Code := -1;
    KillTimer(0, Timer);
    ShowProgress;
    if Code <> 0 then
      MsgBox('The installation did not finish (an internet connection is needed). The messages are in ' +
             ExpandConstant('{app}\install.log') + '. Run this setup again to retry.', mbError, MB_OK);
  end;
end;

function UpdateReadyMemo(Space, NewLine, MemoUserInfoInfo, MemoDirInfo, MemoTypeInfo, MemoComponentsInfo,
  MemoGroupInfo, MemoTasksInfo: String): String;
begin
  Result := MemoDirInfo + NewLine + NewLine + 'Installation steps:' + NewLine + Space + InstallArguments;
  if MemoTasksInfo <> '' then
    Result := Result + NewLine + NewLine + MemoTasksInfo;
end;
