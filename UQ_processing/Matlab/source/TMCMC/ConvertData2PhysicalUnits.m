% function ConvertData2PhysicalUnits

% copy results of inferene here:


% Patient = 1;
%L = 25.6;
% best     = [0.0018833050427   0.0288319149179   273.1305177658181   0.6613240456011 0.6916730293767   0.5517100652865   0.1107012416918   0.9586122727895 0.6925187575929   0.5213808622745   0.0506570811165  -230.8909700000000]
% meanData = [ 0.0033377792556   0.0539495649369   192.6492538009270   0.6610707846963 0.6918143696683   0.5501647681315   0.1093939740974   0.9512963188920  0.7119513264799   0.5156379216779   0.0535420410999  -235.6620764798610]
% stdData  = [ 0.0019020993637   0.0306124161039   107.3959388956695   0.0016811280315  0.0015279281938   0.0016847537845   0.0131130333000   0.0403149448507 0.0293594892376   0.0163268718128   0.0022028024133   2.3904320126141 ]

% 
% Patient = 7;
%L = 25.6;
% best     =[0.0006005514828   0.0238130900919   337.5824782298440   0.4225730668277 0.6179028701972   0.6794223159029   0.1602250047303   0.8457813717993  0.5585570363624   0.4361128708549   0.0507  -340.5082850000000];
% meanData =[0.0009028124893   0.0357381940647   407.1009591298574   0.4228936336877 0.6171250701941   0.6790941193711   0.1608859172769   0.8584818543687  0.5772038409151   0.4474750421409   0.0529322259155  -346.8389858696766];
% stdData  =[0.0008558726794   0.0336298318620   257.5284375113604   0.0014920091925  0.0014341180072   0.0019943717949   0.0141182287961   0.0547641104799  0.0265798883287   0.0163842172645   0.0025248963410   3.1557361532395];

%  Patient = 11;
%  L = 25.6;

 % 4K
% best =     [0.000893100558   0.006449208259   1809.148154216019 0.354802755465   0.556997135993   0.776283236687  0.194106994962   0.742004175122   0.733878243465  0.307727510365   0.050026075219  -2262.476994000000];
% meanData = [0.001112510868   0.007819084917   1571.611005681791 0.349713620255   0.560385703448   0.773509823727 0.192303804270   0.751604428991   0.747495206683  0.309068579697   0.053242138716  -2269.731473935190];
% stdData =  [0.0003436365306   0.0022835216712   389.2764735766627 0.0039004549521   0.0033995989131   0.0024017554460  0.0052388484748   0.0247856130787   0.0171346050184 0.0110484853483   0.0021528650999   3.2905689186828];

% 6K
% best     = [ 0.32566088976   0.0133381534040   822.5525607088342   0.3386667042329   0.5754106061544   0.7730853285253   0.1572646894220   0.6003779211718 0.6224679656862   0.2902659081652   0.0527699065018  -513.9943980000001];
% meanData = [ 0.46399446935   0.0195812540142   772.4756569009106   0.3399382101427   0.5730712535073   0.7774247193704   0.1662724799737   0.6086694566007 0.6298452079974   0.2982349141017   0.0536088464754  -525.0646490825047];
% stdData  = [ 0.42619745356   0.0174208138573   344.7616632677367   0.0042604058862   0.0031586174229   0.0051063908934   0.0072861748283   0.0050930289735 0.0129926745397   0.0088791509967   0.0032987942187   4.1743620626261];


% Patient = 22;
% L = 25.6;
% 4K
% best = [0.011662100667   0.073704201927   138.518274746381   0.646292744466 0.471314537437   0.632733952051   0.178890361701   0.882590415042   0.555135421705   0.346938061176   0.050143544907  -2468.336935000000];
% meanData = [0.004384337191   0.028508547051   546.841537812448   0.646062694492  0.470767691846   0.631749989142   0.184423781328   0.890146054831 0.555309308954   0.348741992844   0.050615198380  -2477.004862748143];
% stdData =[0.0038466622007   0.0254904951384   292.8739876592243   0.0008783554119   0.0007521821558   0.0008095034324   0.0038917035034   0.0142099906045  0.0071064410818   0.0031754552543   0.0004045219883   2.2622865219893];


% 6K 
% best    =  [ 0.001241606644    0.008461051593    1248.948234044119   0.645273133782   0.470134616432    0.624025479777    0.173930168793    0.867798219561   0.601022014115   0.352924042257   0.056654727567  -622.376392000000];
% meanData = [ 0.0038063739500   0.0252166383463   672.9740729883814   0.6447697296450  0.4696619326955   0.6248234564800   0.1720081437778   0.8648809911100  0.6064223062898  0.3613758981144  0.0589232886156  -629.3807400772157];
% stdData =  [ 0.0034315124003   0.0214494156428   436.3092905851151   0.0014014111455  0.0012287004034   0.0013482306759   0.0063696248246   0.0284740547633 0.0049981972232   0.0123597043837  0.0037701086509   3.0914467302580];

%Synthetic small
% L = 21.7;
% 
% best     = [0.0012255479124   0.0295994434680   271.1326786610761   0.3148831561813  0.6701239890387   0.4997306359439   0.0280459261612   0.9884378981863   0.8440068507555   0.4874448058985   0.0501601596713  -125.4801460000000];
% meanData = [0.0015427796284   0.0360958137415   253.8467984931114   0.3149682669743  0.6698286893060   0.4992049853459   0.0294237207053   0.9774877608949  0.8324511927874   0.4840847212395   0.0504914959984  -134.8427573311013 ];
% stdData  = [0.0005295962858   0.0124737411225   107.5892009021741   0.0004681821289  0.0004304265926   0.0007674401193   0.0022116295243   0.0139206983463   0.0131529608469   0.0075686643112   0.0004157419384   3.7938547696748 ];


% Synthetic Big
 L = 25.6;
% 4 K
% best = [0.0016767196308   0.0101110624912   887.0072511444196   0.5995983450818   0.4497626746289   0.6503288134006   0.0346850286873   0.9088966616145  0.7806217530226   0.4128300134212   0.0511600652879  -304.5673060000000];
% meanData = [0.001334099658   0.007858077953   1287.386095146544   0.599609606088  0.449866611030   0.650598266295   0.035536655841   0.897751229267   0.773569241513   0.412613326926   0.050690302231  -323.623445763731];
% stdData = [0.0004633285675   0.0028337695288   454.7084953129391   0.0009307650062   0.0003543936317   0.0006769015759   0.0011594580217   0.0091207591309  0.0060261390752   0.0048236957137   0.0003408135667   5.6341234381207];

% 6 K

best     = [ 0.115310812798   0.006632647972   1345.528115895332   0.599908494098 0.449964754839   0.649784067012   0.022946805255   0.900567132243 0.725133606100   0.429913174476   0.054047750596   673.489420000000];
meanData = [ 0.146864538752   0.008356090511   1129.185983248272   0.599889391902 0.449562701081   0.650003123508   0.023973206075   0.901089185642 0.706055916530   0.434529544679   0.053836194567   654.080304544806];
stdData  = [ 0.0371169877151  0.0021275829371  251.4558429094923   0.0006406075519  0.0003929563403  0.0003838927554  0.0008744060950   0.0102529903098 0.0241179512359  0.0089947714696   0.0015612780172   6.0667202429135];
% conversion:
%1) Diffusion from cm2 / day -> mm^2 day
% best(1) = best(1) * 100;
% meanData(1) = meanData(1) * 100;
% stdData(1) = stdData(1) * 100;

%2) IC, convert from normalised units to real one and then to mm: * L * 10
best(4) = best(4) * L * 10;  % ix
best(5) = best(5) * L * 10;  % iy
best(6) = best(6) * L * 10;  % iz

meanData(4) = meanData(4) * L * 10;  % ix
meanData(5) = meanData(5) * L * 10;  % iy
meanData(6) = meanData(6) * L * 10;  % iz

stdData(4) = stdData(4) * L * 10;  % ix
stdData(5) = stdData(5) * L * 10;  % iy
stdData(6) = stdData(6) * L * 10;  % iz


best
meanData
stdData

