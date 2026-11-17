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
Name: "{userprograms}\Foci Analyzer"; Filename: "{app}\.venv\Scripts\foci-analyzer-window.exe"; WorkingDir: "{app}"; IconFilename: "{app}\program\installer\FociAnalyzer.ico"; Comment: "Foci in nuclei and cells"
Name: "{userdesktop}\Foci Analyzer"; Filename: "{app}\.venv\Scripts\foci-analyzer-window.exe"; WorkingDir: "{app}"; IconFilename: "{app}\program\installer\FociAnalyzer.ico"; Tasks: desktopicon

[UninstallDelete]
; the environments, Python, the download cache and the log were made by install.py, not copied by the setup
Type: filesandordirs; Name: "{app}"

[Run]
Filename: "{app}\.venv\Scripts\foci-analyzer-window.exe"; Description: "Start Foci Analyzer"; Flags: postinstall nowait skipifsilent unchecked

[Code]
const
  NetCount = 5;
var
  NetPage: TInputOptionWizardPage;
  OwnPage: TInputFileWizardPage;
  Names: array[0..4] of String;

procedure InitializeWizard;
var
  i: Integer;
begin
  Names[0] := 'stardist'; Names[1] := 'cellpose'; Names[2] := 'instanseg'; Names[3] := 'spotiflow'; Names[4] := 'trackastra';
  NetPage := CreateInputOptionPage(wpSelectTasks, 'Neural networks',
    'Which networks should be installed?',
    'Foci Analyzer works without them (threshold segmentation, watershed foci), but they find nuclei and foci ' +
    'better. They need about 7 GB together and run on an NVIDIA graphics card when there is one. ' +
    'You can run this setup again later to add them.', False, False);
  NetPage.Add('StarDist: nuclei, fast and robust (recommended)');
  NetPage.Add('Cellpose: nuclei or whole cells, also in 3D');
  NetPage.Add('InstanSeg: nuclei or whole cells, fast; also made for tissue');
  NetPage.Add('Spotiflow: foci detection by a network');
  NetPage.Add('Trackastra: tracking of cells in time-lapses');
  for i := 0 to NetCount - 1 do
    NetPage.Values[i] := True;

  OwnPage := CreateInputFilePage(NetPage.ID, 'Your own environments (optional)',
    'Already have one of the networks installed in a Python environment of your own?',
    'Give the python.exe of that environment (for example C:\PythonProjects\Cellpose4\.venv\Scripts\python.exe). ' +
    'That network is then not installed again, and Foci Analyzer uses yours. Leave a box empty otherwise. ' +
    'This can be changed later on the Settings page of Foci Analyzer.');
  for i := 0 to NetCount - 1 do
    OwnPage.Add(Names[i] + ':', 'python.exe|python.exe|All files|*.*', '.exe');
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
    if NetPage.Values[i] and (Trim(OwnPage.Values[i]) = '') then
    begin
      if Nets <> '' then Nets := Nets + ',';
      Nets := Nets + Names[i];
    end;
  if Nets = '' then Nets := 'none';
  Result := '--networks ' + Nets;
  for i := 0 to NetCount - 1 do
    if Trim(OwnPage.Values[i]) <> '' then
      Result := Result + ' --own ' + Quote(Names[i] + '=' + Trim(OwnPage.Values[i]));
end;

function NextButtonClick(CurPageID: Integer): Boolean;
var
  i: Integer;
begin
  Result := True;
  if CurPageID = OwnPage.ID then
    for i := 0 to NetCount - 1 do
      if (Trim(OwnPage.Values[i]) <> '') and not FileExists(Trim(OwnPage.Values[i])) then
      begin
        MsgBox('There is no file ' + Trim(OwnPage.Values[i]) + ' (for ' + Names[i] + ').', mbError, MB_OK);
        Result := False;
        Exit;
      end;
end;

procedure CurStepChanged(CurStep: TSetupStep);
var
  Code: Integer;
  Script: String;
begin
  if CurStep = ssPostInstall then
  begin
    WizardForm.StatusLabel.Caption := 'Downloading Python, the program''s packages and the networks; ' +
      'this can take a while (the black window shows what happens)...';
    Script := ExpandConstant('{app}\program\installer\run-install.cmd');
    if not Exec(ExpandConstant('{cmd}'), '/c "' + Quote(Script) + ' ' + InstallArguments + '"', '', SW_SHOW,
                ewWaitUntilTerminated, Code) or (Code <> 0) then
      MsgBox('The installation did not finish (an internet connection is needed). The log is in ' +
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
