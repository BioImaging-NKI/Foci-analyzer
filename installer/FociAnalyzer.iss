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
DisableDirPage=no
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
; Inno Setup 6.7+ enables Windows' RedirectionGuard, and the programs Setup starts inherit it. uv links
; python\cpython-3.12-... to the patch version with a junction made without administrator rights, which
; RedirectionGuard refuses to follow: "untrusted mount point (os error 448)". So it stays off.
RedirectionGuard=no

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
Filename: "{app}\Foci Analyzer.exe"; Parameters: "--after-update"; Description: "Start Foci Analyzer"; Flags: postinstall nowait skipifsilent unchecked

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
  { an installation that is there already (from the Apps entry the setup wrote) }
  PrevDir, PrevVersion: String;
  ExistPage: TInputOptionWizardPage;
  Installed: array[0..5] of Boolean;
  { the networks the user took from an environment of their own (program.json of the window) }
  OwnBefore: array[0..5] of String;
  Removing: Boolean;

{ While the installation runs, Setup waits in its own loop and handles the window's messages itself: a script
  Exec that waits disables the whole wizard (it cannot be moved, scrolled or clicked). Setup is 32-bit. }
type
  TWinMsg = record
    hwnd: Longword;
    message: Longword;
    wParam: Longint;
    lParam: Longint;
    time: Longword;
    ptX: Longint;
    ptY: Longint;
  end;

function PeekMessage(var Msg: TWinMsg; hWnd, MsgFilterMin, MsgFilterMax, RemoveMsg: Longword): Bool;
external 'PeekMessageW@user32.dll stdcall';
function TranslateMessage(var Msg: TWinMsg): Bool;
external 'TranslateMessage@user32.dll stdcall';
function DispatchMessage(var Msg: TWinMsg): Longint;
external 'DispatchMessageW@user32.dll stdcall';

procedure HandleMessages;
var
  Msg: TWinMsg;
begin
  while PeekMessage(Msg, 0, 0, 0, 1) do  { PM_REMOVE }
  begin
    TranslateMessage(Msg);
    DispatchMessage(Msg);
  end;
end;

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
  { only what is not installed yet counts }
  Total := 0;
  if PrevDir = '' then Total := SizeProgram;
  Torch := False;
  for i := 0 to NetCount - 1 do
    if NetPage.Values[i] and not Installed[i] then
    begin
      Total := Total + Sizes[i];
      if (i >= 1) and (i <= 4) then Torch := True;
      if i = 5 then Total := Total + TorchSize;
    end;
  if Torch and not (Installed[1] or Installed[2] or Installed[3] or Installed[4]) then Total := Total + TorchSize;
  if PrevDir = '' then
    TotalLabel.Caption := 'Together with the program: about ' + GB(Total) + ' on disk.'
  else
    TotalLabel.Caption := 'Added to the installation: about ' + GB(Total) + ' on disk.';
  if NetPage.Values[5] then
    TotalLabel.Caption := TotalLabel.Caption + ' Cellpose 4 downloads its model (CellposeSAM v2, 1.2 GB) when it ' +
      'is first used.';
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

function HasPackage(Env, Dist: String): Boolean;
var
  R: TFindRec;
begin
  Result := FindFirst(PrevDir + '\' + Env + '\Lib\site-packages\' + Dist + '-*.dist-info', R);
  if Result then FindClose(R);
end;

{ the installation this setup made earlier: where, which version, which networks }
procedure FindPrevious;
var
  Key: String;
begin
  Key := 'Software\Microsoft\Windows\CurrentVersion\Uninstall\{6C1F7E52-3A0B-4C55-9E4D-F0C1A11A1A11}_is1';
  PrevDir := '';
  PrevVersion := '';
  if RegQueryStringValue(HKCU, Key, 'InstallLocation', PrevDir) then
  begin
    PrevDir := RemoveBackslashUnlessRoot(PrevDir);
    if not FileExists(PrevDir + '\unins000.exe') then PrevDir := '';
  end;
  if PrevDir = '' then Exit;
  RegQueryStringValue(HKCU, Key, 'DisplayVersion', PrevVersion);
  Installed[0] := HasPackage('stardist-env', 'stardist');
  Installed[1] := HasPackage('models-env', 'cellpose');
  Installed[2] := HasPackage('models-env', 'instanseg_torch');
  Installed[3] := HasPackage('models-env', 'spotiflow');
  Installed[4] := HasPackage('models-env', 'trackastra');
  Installed[5] := HasPackage('cellpose4-env', 'cellpose');
end;

{ "Name": "C:\\path\\python.exe" in program.json (written by Python's json module) }
function JsonString(Text, Name: String): String;
var
  P: Integer;
begin
  Result := '';
  P := Pos('"' + Name + '":', Text);
  if P = 0 then Exit;
  Text := Copy(Text, P + Length(Name) + 3, Length(Text));
  Text := TrimLeft(Text);
  if (Text = '') or (Text[1] <> '"') then Exit;
  P := 2;
  while (P <= Length(Text)) and (Text[P] <> '"') do
  begin
    if Text[P] = '\' then
    begin
      P := P + 1;
      if P > Length(Text) then Exit;
    end;
    Result := Result + Text[P];
    P := P + 1;
  end;
end;

procedure FindOwn;
var
  Text: AnsiString;
  i, P: Integer;
begin
  if not LoadStringFromFile(ExpandConstant('{%USERPROFILE}\.foci_analyzer\program.json'), Text) then Exit;
  P := Pos('"pythons":', Text);
  if P = 0 then Exit;
  Text := Copy(Text, P, Length(Text));
  P := Pos('}', Text);
  if P > 0 then Text := Copy(Text, 1, P);
  for i := 0 to NetCount - 1 do
    OwnBefore[i] := JsonString(Text, Names[i]);
end;

procedure InitializeWizard;
var
  i: Integer;
  Torch, Which: String;
begin
  FindPrevious;
  if PrevDir <> '' then
  begin
    if PrevVersion = '{#AppVersion}' then
      Which := 'This version (' + PrevVersion + ') is installed'
    else
      Which := 'Foci Analyzer ' + PrevVersion + ' is installed';
    ExistPage := CreateInputOptionPage(wpWelcome, 'Foci Analyzer is installed already',
      Which + ' in ' + PrevDir + '.',
      'This setup has version {#AppVersion}. Updating keeps your settings, experiments, results and the networks ' +
      'that are installed; on the next pages you can add networks. Removing takes Foci Analyzer and its networks ' +
      'off this computer (your settings and results stay).', True, False);
    ExistPage.Add('Update it, or add networks');
    ExistPage.Add('Remove it');
    ExistPage.SelectedValueIndex := 0;
  end;
  Names[0] := 'stardist'; Names[1] := 'cellpose'; Names[2] := 'instanseg'; Names[3] := 'spotiflow'; Names[4] := 'trackastra'; Names[5] := 'cellpose4';
  Shown[0] := 'StarDist'; Shown[1] := 'Cellpose 3'; Shown[2] := 'InstanSeg'; Shown[3] := 'Spotiflow'; Shown[4] := 'Trackastra'; Shown[5] := 'Cellpose 4';
  FindOwn;
  Sizes[0] := 16; Sizes[1] := 1; Sizes[2] := 3; Sizes[3] := 2; Sizes[4] := 2; Sizes[5] := 1;
  HasGpu := FileExists(ExpandConstant('{sys}\nvcuda.dll')) or FileExists(ExpandConstant('{sys}\nvidia-smi.exe'));
  if HasGpu then
    Torch := 'PyTorch for the NVIDIA graphics card (' + GB(TorchSize) + ')'
  else
    Torch := 'PyTorch for the processor (' + GB(TorchSize) + '; no NVIDIA graphics card found)';

  NetPage := CreateInputOptionPage(wpSelectTasks, 'Neural networks',
    'Which networks should be installed?',
    'Foci Analyzer works without them (threshold segmentation, watershed foci), but they find nuclei and foci ' +
    'better. Cellpose 3, InstanSeg, Spotiflow and Trackastra share one ' + Torch + '; Cellpose 4 has its own. ' +
    'StarDist runs on the processor. Networks you leave unticked can be taken from a Python environment of your ' +
    'own on the next page, or added later.', False, False);
  NetPage.Add('StarDist: nuclei, fast and robust (recommended) - ' + GB(Sizes[0]));
  NetPage.Add('Cellpose 3: nuclei or whole cells, also in 3D - ' + GB(Sizes[1]) + ' + PyTorch');
  NetPage.Add('InstanSeg: nuclei or whole cells, fast; also made for tissue - ' + GB(Sizes[2]) + ' + PyTorch');
  NetPage.Add('Spotiflow: foci detection by a network - ' + GB(Sizes[3]) + ' + PyTorch');
  NetPage.Add('Trackastra: tracking of cells in time-lapses - ' + GB(Sizes[4]) + ' + PyTorch');
  NetPage.Add('Cellpose 4: CellposeSAM, the newest model; slow without a modern graphics card - ' + GB(Sizes[5] + TorchSize));
  for i := 0 to NetCount - 2 do
    NetPage.Values[i] := PrevDir = '';
  for i := 0 to NetCount - 1 do
    if Installed[i] then
    begin
      NetPage.Values[i] := True;
      NetPage.CheckListBox.ItemEnabled[i] := False;
      NetPage.CheckListBox.ItemCaption[i] := NetPage.CheckListBox.ItemCaption[i] + ' (installed)';
    end;
  { a network from the user's own environment is not ticked; ticking it installs it and uses that instead }
  for i := 0 to NetCount - 1 do
    if OwnBefore[i] <> '' then
    begin
      if not Installed[i] then NetPage.Values[i] := False;
      NetPage.CheckListBox.ItemCaption[i] := NetPage.CheckListBox.ItemCaption[i] + ' - now from your own environment';
    end;
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
  OwnIntro.Height := ScaleY(100);
  OwnIntro.Anchors := [akLeft, akTop, akRight];
  OwnIntro.Caption := 'These networks are not installed. If one of them is installed in a Python environment of ' +
    'your own, give the python.exe of that environment: Foci Analyzer then uses it. A venv (for example ' +
    'C:\PythonProjects\Cellpose4\.venv\Scripts\python.exe) and a conda environment (for example ' +
    'C:\Users\you\miniconda3\envs\cellpose\python.exe) both work. Leave a box empty otherwise; this can be ' +
    'changed later on the Settings page.';
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
    OwnEdit[i].Text := OwnBefore[i];
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
  if (PageID = wpSelectDir) and (PrevDir <> '') then
    Result := True;  { an update goes into the folder that is there }
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
      Result := Result + ' --own ' + Quote(Names[i] + '=' + OwnPython(i))
    else if (OwnBefore[i] <> '') and not Installed[i] then
      Result := Result + ' --own ' + Names[i] + '=';  { no longer the user's own: the installed one }
end;

function NextButtonClick(CurPageID: Integer): Boolean;
var
  i: Integer;
begin
  Result := True;
  if (PrevDir <> '') and (CurPageID = ExistPage.ID) and (ExistPage.SelectedValueIndex = 1) then
  begin
    { the uninstaller of that installation; this setup then closes }
    if Exec(PrevDir + '\unins000.exe', '', '', SW_SHOW, ewNoWait, i) then
    begin
      Removing := True;
      WizardForm.Close;
    end
    else
      MsgBox('The uninstaller could not be started: ' + SysErrorMessage(i), mbError, MB_OK);
    Result := False;
    Exit;
  end;
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

{ the last lines of the installation: from the window, or, when Python did not even start, from the output file }
function LastMessages: String;
var
  Lines: TArrayOfString;
  i, First: Integer;
begin
  Result := '';
  if LogMemo.Lines.Count > 0 then
  begin
    First := LogMemo.Lines.Count - 8;
    if First < 0 then First := 0;
    for i := First to LogMemo.Lines.Count - 1 do
      Result := Result + LogMemo.Lines[i] + #13#10;
  end
  else if LoadStringsFromFile(ExpandConstant('{app}\install-output.log'), Lines) then
  begin
    First := GetArrayLength(Lines) - 8;
    if First < 0 then First := 0;
    for i := First to GetArrayLength(Lines) - 1 do
      Result := Result + Lines[i] + #13#10;
  end;
  if Trim(Result) = '' then Result := '(none)';
end;

{ A Foci Analyzer that runs from the installation folder (the program, or its network processes) keeps its
  files in use, so they cannot be replaced or removed: it is closed first. }
function RunningCount(Dir: String; Stop: Boolean): Integer;
var
  Command: String;
begin
  Command := '$d = ''' + RemoveBackslashUnlessRoot(Dir) + '\''; ' +
    '$p = @(Get-CimInstance Win32_Process | Where-Object { $_.ExecutablePath -and ' +
    '$_.ExecutablePath.StartsWith($d, [StringComparison]::OrdinalIgnoreCase) }); ';
  if Stop then
    Command := Command + '$p | ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }; ' +
      'Start-Sleep -Seconds 2; ';
  Command := Command + 'exit $p.Count';
  if not Exec('powershell.exe', '-NoProfile -NonInteractive -ExecutionPolicy Bypass -Command "' + Command + '"', '',
              SW_HIDE, ewWaitUntilTerminated, Result) then
    Result := 0;
end;

function CloseRunning(Dir: String; Silent: Boolean): Boolean;
begin
  Result := True;
  if not DirExists(Dir) or (RunningCount(Dir, False) = 0) then Exit;
  if not Silent then
    Result := MsgBox('Foci Analyzer is running from ' + Dir + '. Setup closes it first (your experiment and ' +
                     'settings are kept; its browser tab waits and comes back when Foci Analyzer starts again).', mbConfirmation, MB_OKCANCEL) = IDOK;
  if Result then
    RunningCount(Dir, True);
end;

function PrepareToInstall(var NeedsRestart: Boolean): String;
begin
  Result := '';
  if not CloseRunning(WizardDirValue, WizardSilent) then
    Result := 'Foci Analyzer is still running. Close it and run this setup again.';
end;

function InitializeUninstall: Boolean;
begin
  Result := CloseRunning(ExpandConstant('{app}'), UninstallSilent);
end;

procedure CurStepChanged(CurStep: TSetupStep);
var
  Code, Ticks: Integer;
  Script, ExitFile, Text: AnsiString;
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
    { run-install.cmd writes its exit code to install-exit.txt when it is done }
    ExitFile := ExpandConstant('{app}\install-exit.txt');
    DeleteFile(ExitFile);
    Code := -1;
    if Exec(ExpandConstant('{cmd}'), '/c "' + Quote(Script) + ' ' + InstallArguments + ' --progress ' +
            Quote(ProgressFile) + '"', '', SW_HIDE, ewNoWait, Code) then
    begin
      Ticks := 0;
      while not FileExists(ExitFile) do
      begin
        HandleMessages;
        Sleep(40);
        Ticks := Ticks + 1;
        if Ticks mod 12 = 0 then ShowProgress;
        { about once a minute: has the installation stopped without saying so? (nothing of it runs any more) }
        if (Ticks mod 1500 = 0) and not FileExists(ExitFile) and (RunningCount(ExpandConstant('{app}'), False) = 0) then
        begin
          Sleep(2000);
          if not FileExists(ExitFile) then Break;
        end;
      end;
      Sleep(200);
      if LoadStringFromFile(ExitFile, Text) then
        Code := StrToIntDef(Trim(String(Text)), -1);
    end;
    ShowProgress;
    if Code <> 0 then
      MsgBox('The installation did not finish. The last messages:' + #13#10#13#10 + LastMessages + #13#10#13#10 +
             'All messages are in ' + ExpandConstant('{app}\install-output.log') + '. Run this setup again to retry.',
             mbError, MB_OK);
  end;
end;

function UpdateReadyMemo(Space, NewLine, MemoUserInfoInfo, MemoDirInfo, MemoTypeInfo, MemoComponentsInfo,
  MemoGroupInfo, MemoTasksInfo: String): String;
begin
  Result := MemoDirInfo + NewLine + NewLine + 'Installation steps:' + NewLine + Space + InstallArguments;
  if MemoTasksInfo <> '' then
    Result := Result + NewLine + NewLine + MemoTasksInfo;
end;

procedure CancelButtonClick(CurPageID: Integer; var Cancel, Confirm: Boolean);
begin
  if Removing then Confirm := False;  { closed because the uninstaller took over }
end;
