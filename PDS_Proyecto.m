%% PROYECTO FINAL - LIMPIEZA DE VOZ CON FILTROS FIR E IIR
clear; close all; clc;
archivoA = 'PDS_voz_60Hz_4s.wav';
archivoB = 'PDS_voz_fuente_simulada_4s.wav';
ESCUCHAR = true;                 % Escucha original, FIR e IIR finales.
ESCUCHAR_CANDIDATOS_B = false;   % Activar para comparar cortes de 2--5 kHz.
tramo = [3.8 3.9];               % Tramo de baja actividad para medir RMS.

%% 1. CARGA, REVISION Y GRAFICAS TEMPORALES
[xA,fsA] = audioread(archivoA); xA = mono(xA);
[xB,fsB] = audioread(archivoB); xB = mono(xB);
assert(fsA == fsB,'Los audios deben tener la misma frecuencia de muestreo.');
fs = fsA;
assert(fs == 16000,'Se esperaban audios muestreados a 16 kHz.');

tA = (0:numel(xA)-1)'/fs; tB = (0:numel(xB)-1)'/fs;
idxA = tA >= tramo(1) & tA < tramo(2);
idxB = tB >= tramo(1) & tB < tramo(2);

Audio = {'A: interferencia de 60 Hz';'B: fuente de agua simulada'};
Muestras = [numel(xA);numel(xB)];
Duracion_s = Muestras/fs;
Pico = [max(abs(xA));max(abs(xB))];
disp(table(Audio,repmat(fs,2,1),Muestras,Duracion_s,Pico, ...
    'VariableNames',{'Audio','Fs_Hz','Muestras','Duracion_s','Pico'}));

figure('Name','Diagnostico temporal','Color','w');
tiledlayout(2,2,'TileSpacing','compact');
nexttile; plot(tA,xA); grid on; title('Audio A completo');
xlabel('Tiempo [s]'); ylabel('Amplitud');
nexttile; plot(tA(idxA),xA(idxA)); grid on; title('Audio A: 3.8--3.9 s');
xlabel('Tiempo [s]'); ylabel('Amplitud');
nexttile; plot(tB,xB); grid on; title('Audio B completo');
xlabel('Tiempo [s]'); ylabel('Amplitud');
nexttile; plot(tB(idxB),xB(idxB)); grid on; title('Audio B: 3.8--3.9 s');
xlabel('Tiempo [s]'); ylabel('Amplitud');

%% 2. DIAGNOSTICO FFT POR SEGMENTOS CON VENTANA HANN
% Hann aqui reduce fuga espectral antes de cada FFT. No es todavia la
% ventana que truncara las respuestas impulsionales ideales de los FIR.
longitudes = [1024 4096];
Duracion_ms = 1000*longitudes'/fs;
Resolucion_Hz = fs./longitudes';
disp('COMPROMISO TIEMPO--FRECUENCIA');
disp(table(longitudes',Duracion_ms,Resolucion_Hz, ...
    'VariableNames',{'N','Duracion_ms','Delta_f_Hz'}));

for k = 1:numel(longitudes)
    L = longitudes(k);
    [fL,PA{k},PAseg{k}] = psdSegmentada(xA,fs,L); %#ok<SAGROW>
    [~, PB{k},PBseg{k}] = psdSegmentada(xB,fs,L); %#ok<SAGROW>
end

% Espectros por segmento y promedio: vista completa y bajas frecuencias.
k = numel(longitudes); f = fL; nMostrar = 6;
figure('Name','FFT por segmentos y promedio','Color','w');
tiledlayout(2,2,'TileSpacing','compact');
graficaSegmentos(f,PAseg{k},PA{k},[0 fs/2],nMostrar,'Audio A: espectro completo');
graficaSegmentos(f,PAseg{k},PA{k},[0 300], nMostrar,'Audio A: bajas frecuencias');
graficaSegmentos(f,PBseg{k},PB{k},[0 fs/2],nMostrar,'Audio B: espectro completo');
graficaSegmentos(f,PBseg{k},PB{k},[0 300], nMostrar,'Audio B: bajas frecuencias');

% Comparacion de longitudes y de ventana rectangular contra Hann.
[fR,AR] = periodogramaTramo(xA,fs,4096,'rectangular');
[~, AH] = periodogramaTramo(xA,fs,4096,'hann');
[~, BR] = periodogramaTramo(xB,fs,4096,'rectangular');
[~, BH] = periodogramaTramo(xB,fs,4096,'hann');

figure('Name','Longitud y fuga espectral','Color','w');
tiledlayout(2,2,'TileSpacing','compact');
nexttile; hold on;
for k=1:numel(longitudes), plot((0:numel(PA{k})-1)'*fs/longitudes(k),dB10(PA{k})); end
xlim([0 300]); grid on; title('Audio A: longitud de segmento');
legend(compose('N=%d',longitudes),'Location','best'); xlabel('Hz'); ylabel('PSD [dB/Hz]');
nexttile; hold on;
for k=1:numel(longitudes), plot((0:numel(PB{k})-1)'*fs/longitudes(k),dB10(PB{k})); end
xlim([0 fs/2]); grid on; title('Audio B: longitud de segmento');
legend(compose('N=%d',longitudes),'Location','best'); xlabel('Hz'); ylabel('PSD [dB/Hz]');
nexttile; plot(fR,dB10(AR),fR,dB10(AH)); xlim([0 300]); grid on;
title('Audio A: rectangular vs Hann'); legend('Rectangular','Hann'); xlabel('Hz'); ylabel('PSD [dB/Hz]');
nexttile; plot(fR,dB10(BR),fR,dB10(BH)); xlim([0 fs/2]); grid on;
title('Audio B: rectangular vs Hann'); legend('Rectangular','Hann'); xlabel('Hz'); ylabel('PSD [dB/Hz]');

[~,picoA] = max(PA{end}(f <= 100)); fBajo = f(f <= 100);
fprintf('\nDIAGNOSTICO\nPico dominante del Audio A: %.2f Hz\n',fBajo(picoA));
fprintf(['Audio A: interferencia tonal estrecha -> rechazo 50--70 Hz.\n' ...
         'Audio B: ruido de banda ancha solapado con la voz -> reduccion parcial.\n']);

%% 3. AUDIO A - PRUEBA Y DISENO FIR RECHAZA-BANDA
bandaA = [50 70]; ordenesA = [1000 2000 4000];
at60 = zeros(size(ordenesA)); retardoA_ms = at60; HA_cand = cell(size(ordenesA));
for k=1:numel(ordenesA)
    h = firHann('stop',fs,bandaA,ordenesA(k));
    [H,fH] = freqz(h,1,65536,fs); HA_cand{k}=dB20(H);
    at60(k)=interp1(fH,HA_cand{k},60); retardoA_ms(k)=500*ordenesA(k)/fs;
end
disp('CANDIDATOS FIR DEL AUDIO A');
disp(table(ordenesA',ordenesA'+1,at60',retardoA_ms', ...
    'VariableNames',{'Orden','Coeficientes','Respuesta_60Hz_dB','Retardo_ms'}));

figure('Name','Ordenes FIR del Audio A','Color','w'); hold on;
for k=1:numel(ordenesA), plot(fH,HA_cand{k},'LineWidth',1.1); end
xline(50,':k'); xline(60,'--r'); xline(70,':k'); xlim([20 110]); ylim([-100 5]); grid on;
legend(compose('Orden %d',ordenesA),'Location','best');
xlabel('Frecuencia [Hz]'); ylabel('Magnitud [dB]'); title('Selección del orden del FIR A');

% Seleccion: orden 4000. Fue el candidato que hizo inaudible el tono en
% las pausas; su costo es 4001 multiplicaciones/muestra y 125 ms causales.
ordenA = 4000; hA = firHann('stop',fs,bandaA,ordenA); DA = ordenA/2;
yA_FIR_causal = filter(hA,1,xA);
yA_FIR = firCompensado(hA,xA);
[HA,fHA] = freqz(hA,1,65536,fs);
[gdA,fgA] = grpdelay(hA,1,4096,fs);

fprintf('\nFIR A\nH_A(z)=sum_{k=0}^{%d} h_A[k] z^{-k}\n',ordenA);
fprintf('y[n]=sum_{k=0}^{%d} h_A[k]x[n-k]\n',ordenA);
fprintf('Fase lineal; retardo causal = %d muestras = %.1f ms.\n',DA,1000*DA/fs);

figure('Name','Diseno FIR del Audio A','Color','w');
tiledlayout(2,2,'TileSpacing','compact');
nexttile; plot(0:ordenA,hA); grid on; title('h_A[n]'); xlabel('n'); ylabel('Amplitud');
nexttile; plot(fHA,dB20(HA)); xlim([20 110]); ylim([-100 5]); grid on;
xline(60,'--r'); title('Magnitud'); xlabel('Hz'); ylabel('dB');
nexttile; plot(fHA,unwrap(angle(HA))); xlim([0 110]); grid on;
title('Fase lineal'); xlabel('Hz'); ylabel('rad');
nexttile; plot(fgA,gdA); xlim([0 500]); ylim([DA-5 DA+5]); grid on;
title('Retardo de grupo'); xlabel('Hz'); ylabel('muestras');

%% 4. AUDIO A - IIR NOTCH, PLANO Z Y COMPARACION
Q = 30; w0 = 60/(fs/2); [bAIIR,aAIIR] = iirnotch(w0,w0/Q);
yA_IIR_causal = filter(bAIIR,aAIIR,xA);
yA_IIR = filtfilt(bAIIR,aAIIR,xA); % Resultado fuera de linea, fase cero.
[HAIIR,fAIIR] = freqz(bAIIR,aAIIR,65536,fs);
[gdAIIR,fgAIIR] = grpdelay(bAIIR,aAIIR,4096,fs);
pAIIR = roots(aAIIR); radioA = max(abs(pAIIR));

fprintf('\nIIR A\nH(z)=(%.8f %+.8fz^-1 %+.8fz^-2)/(1 %+.8fz^-1 %+.8fz^-2)\n', ...
    bAIIR(1),bAIIR(2),bAIIR(3),aAIIR(2),aAIIR(3));
fprintf(['y[n]=%.8fx[n] %+.8fx[n-1] %+.8fx[n-2] ' ...
         '%+.8fy[n-1] %+.8fy[n-2]\n'], ...
    bAIIR(1),bAIIR(2),bAIIR(3),-aAIIR(2),-aAIIR(3));
fprintf('Radio maximo de polos = %.8f -> estable. filtfilt: fase cero, orden efectivo 4.\n',radioA);

% En el FIR de orden 4000 se muestran los ceros críticos cercanos a 60 Hz;
% calcular las 4000 raíces no añade claridad y es numéricamente costoso.
fCerosA = cerosCriticosFIR(hA,fs,[20 110]);
figure('Name','Plano Z del Audio A','Color','w');
tiledlayout(1,2,'TileSpacing','compact');
nexttile; planoZ(exp(1j*2*pi*fCerosA/fs),[],[0.9994 1.0001 -0.03 0.03]);
title('FIR A: ceros críticos');
nexttile; planoZ(roots(bAIIR),pAIIR,[-1.1 1.1 -1.1 1.1]); title('IIR notch A');

figure('Name','Comparacion FIR e IIR - Audio A','Color','w');
tiledlayout(3,1,'TileSpacing','compact');
nexttile; plot(fHA,dB20(HA),fAIIR,2*dB20(HAIIR),'--');
xlim([20 110]); ylim([-100 5]); grid on; xline(60,':k');
legend('FIR Hann','IIR + filtfilt'); title('Respuesta efectiva'); ylabel('dB');
nexttile; plot(tA(idxA),xA(idxA),'Color',[.65 .65 .65]); hold on;
plot(tA(idxA),yA_FIR(idxA)); plot(tA(idxA),yA_IIR(idxA),'--'); grid on;
legend('Original','FIR','IIR'); title('Comparación temporal'); ylabel('Amplitud');
nexttile; [f0,S0]=psdSegmentada(xA,fs,4096); [~,SF]=psdSegmentada(yA_FIR,fs,4096);
[~,SI]=psdSegmentada(yA_IIR,fs,4096); plot(f0,dB10(S0),f0,dB10(SF),f0,dB10(SI),'--');
xlim([0 300]); grid on; legend('Original','FIR','IIR'); xlabel('Hz'); ylabel('PSD [dB/Hz]');

metA_FIR = metricasA(xA,yA_FIR,fs,idxA);
metA_IIR = metricasA(xA,yA_IIR,fs,idxA);
Metodo = {'FIR Hann';'IIR notch + filtfilt'};
Orden_causal = [ordenA;2]; Orden_efectivo = [ordenA;4];
Atenuacion_60Hz_dB = [metA_FIR(1);metA_IIR(1)];
Cambio_RMS_dB = [metA_FIR(2);metA_IIR(2)];
Correlacion = [metA_FIR(3);metA_IIR(3)];
Mult_causal = [numel(hA);numel(bAIIR)+numel(aAIIR)-1];
disp('COMPARACION FINAL DEL AUDIO A');
disp(table(Metodo,Orden_causal,Orden_efectivo,Atenuacion_60Hz_dB, ...
    Cambio_RMS_dB,Correlacion,Mult_causal));

%% 5. AUDIO B - PRUEBA DE CORTES Y ORDENES FIR
% El Audio A filtrado se usa sólo como referencia aproximada de la misma
% voz. Las correlaciones siguientes NO son una SNR verdadera.
referencia = yA_FIR(1:min(numel(yA_FIR),numel(xB)));
cortesB = [2000 3000 4000 5000]; ordenPrueba = 256;
corrCorte = zeros(size(cortesB)); rmsCorte = corrCorte; hCorte=cell(size(cortesB)); yCorte=hCorte;
figure('Name','Cortes FIR del Audio B','Color','w'); hold on;
for k=1:numel(cortesB)
    hCorte{k}=firHann('low',fs,cortesB(k),ordenPrueba);
    yCorte{k}=firCompensado(hCorte{k},xB);
    [HC,fC]=freqz(hCorte{k},1,32768,fs); plot(fC,dB20(HC),'LineWidth',1.0);
    corrCorte(k)=correlacion(referencia,yCorte{k});
    rmsCorte(k)=cambioRMS(xB,yCorte{k},idxB);
end
xlim([0 fs/2]); ylim([-100 5]); grid on; legend(compose('%g kHz',cortesB/1000));
xlabel('Hz'); ylabel('dB'); title('Frecuencias de corte candidatas');
disp('CANDIDATOS DE CORTE DEL AUDIO B (referencia aproximada, no SNR)');
disp(table(cortesB',repmat(ordenPrueba,numel(cortesB),1),rmsCorte',corrCorte', ...
    'VariableNames',{'Corte_Hz','Orden','Cambio_RMS_dB','Correlacion_auxiliar'}));

if ESCUCHAR_CANDIDATOS_B
    secuencia = xB;
    for k=1:numel(yCorte), secuencia=[secuencia;zeros(round(.35*fs),1);yCorte{k}]; end %#ok<AGROW>
    sound(normaliza(secuencia),fs); pause(numel(secuencia)/fs+.2);
end

% La escucha eligió 2 kHz: fue el audio más limpio, aunque pierde brillo.
fcB = 2000; ordenesB = [128 256 512];
trans40 = zeros(size(ordenesB)); corrOrden=trans40; rmsOrden=trans40; HBcand=cell(size(ordenesB));
figure('Name','Ordenes FIR del Audio B','Color','w'); hold on;
for k=1:numel(ordenesB)
    h = firHann('low',fs,fcB,ordenesB(k)); y = firCompensado(h,xB);
    [H,fH]=freqz(h,1,65536,fs); M=dB20(H); HBcand{k}=M; plot(fH,M,'LineWidth',1.0);
    q=find(fH>=fcB & M<=-40,1); if isempty(q), trans40(k)=NaN; else, trans40(k)=fH(q); end
    corrOrden(k)=correlacion(referencia,y); rmsOrden(k)=cambioRMS(xB,y,idxB);
end
xlim([1500 3000]); ylim([-100 5]); grid on; xline(fcB,'--k');
legend(compose('Orden %d',ordenesB)); xlabel('Hz'); ylabel('dB'); title('Órdenes candidatos');
disp('CANDIDATOS DE ORDEN DEL AUDIO B');
disp(table(ordenesB',ordenesB'+1,500*ordenesB'/fs,trans40',rmsOrden',corrOrden', ...
    'VariableNames',{'Orden','Coeficientes','Retardo_ms','Frecuencia_menos40dB_Hz', ...
    'Cambio_RMS_dB','Correlacion_auxiliar'}));

% Se elige orden 128: escucha equivalente, menor retardo y menor costo.
ordenB = 128; hB = firHann('low',fs,fcB,ordenB); DB = ordenB/2;
yB_FIR_causal = filter(hB,1,xB); yB_FIR = firCompensado(hB,xB);
[HB,fHB]=freqz(hB,1,65536,fs); [gdB,fgB]=grpdelay(hB,1,4096,fs);

fprintf('\nFIR B\nH_B(z)=sum_{k=0}^{%d} h_B[k] z^{-k}\n',ordenB);
fprintf('y[n]=sum_{k=0}^{%d} h_B[k]x[n-k]\n',ordenB);
fprintf('Fase lineal; retardo causal = %d muestras = %.1f ms.\n',DB,1000*DB/fs);

figure('Name','Diseno FIR del Audio B','Color','w');
tiledlayout(2,2,'TileSpacing','compact');
nexttile; stem(0:ordenB,hB,'.'); grid on; title('h_B[n]'); xlabel('n');
nexttile; plot(fHB,dB20(HB)); xlim([0 fs/2]); ylim([-100 5]); grid on;
xline(fcB,'--r'); title('Magnitud'); xlabel('Hz'); ylabel('dB');
nexttile; plot(fHB,unwrap(angle(HB))); xlim([0 fcB]); grid on;
title('Fase lineal'); xlabel('Hz'); ylabel('rad');
nexttile; planoZ(roots(recortaCeros(hB)),[],[-1.2 1.2 -1.2 1.2]); title('Ceros FIR B');

%% 6. AUDIO B - IIR BUTTERWORTH Y COMPARACION FINAL
[bBIIR,aBIIR] = butter(6,fcB/(fs/2),'low');
yB_IIR_causal = filter(bBIIR,aBIIR,xB);
yB_IIR = filtfilt(bBIIR,aBIIR,xB); % Fuera de línea, fase cero.
[HBIIR,fBIIR]=freqz(bBIIR,aBIIR,65536,fs);
[gdBIIR,fgBIIR]=grpdelay(bBIIR,aBIIR,4096,fs);
pBIIR=roots(aBIIR); radioB=max(abs(pBIIR));
fprintf('\nIIR B: Butterworth orden 6, corte 2 kHz. Radio maximo = %.6f -> estable.\n',radioB);
fprintf('Causal: fase no lineal. WAV final: filtfilt, fase cero y orden efectivo 12.\n');

figure('Name','IIR y comparacion final - Audio B','Color','w');
tiledlayout(3,2,'TileSpacing','compact');
nexttile; planoZ(roots(bBIIR),pBIIR,[-1.2 1.2 -1.2 1.2]); title('Butterworth: polos y ceros');
nexttile; plot(fHB,dB20(HB),fBIIR,2*dB20(HBIIR),'--');
xlim([0 fs/2]); ylim([-100 5]); grid on; xline(fcB,':k');
legend('FIR','IIR + filtfilt'); title('Magnitud efectiva'); ylabel('dB');
nexttile; plot(fHB,unwrap(angle(HB)),fBIIR,unwrap(angle(HBIIR)),'--');
xlim([0 4000]); grid on; legend('FIR','IIR causal'); title('Fase causal'); ylabel('rad');
nexttile; plot(fgB,gdB,fgBIIR,gdBIIR,'--'); xlim([0 4000]); grid on;
legend('FIR','IIR causal'); title('Retardo de grupo'); ylabel('muestras');
nexttile; plot(tB(idxB),xB(idxB),'Color',[.7 .7 .7]); hold on;
plot(tB(idxB),yB_FIR(idxB)); plot(tB(idxB),yB_IIR(idxB),'--'); grid on;
legend('Original','FIR','IIR'); title('Comparación temporal'); xlabel('s');
nexttile; [f0,S0]=psdSegmentada(xB,fs,4096); [~,SF]=psdSegmentada(yB_FIR,fs,4096);
[~,SI]=psdSegmentada(yB_IIR,fs,4096); plot(f0,dB10(S0),f0,dB10(SF),f0,dB10(SI),'--');
xlim([0 fs/2]); grid on; legend('Original','FIR','IIR'); title('Comparación espectral'); xlabel('Hz');

metB_FIR = metricasB(xB,yB_FIR,referencia,fs,idxB,fcB);
metB_IIR = metricasB(xB,yB_IIR,referencia,fs,idxB,fcB);
Metodo = {'FIR Hann';'IIR Butterworth + filtfilt'};
Orden_causal = [ordenB;6]; Orden_efectivo = [ordenB;12];
Cambio_RMS_dB = [metB_FIR(1);metB_IIR(1)];
Reduccion_sobre_2kHz_dB = [metB_FIR(2);metB_IIR(2)];
Correlacion_auxiliar = [metB_FIR(3);metB_IIR(3)];
Mult_causal = [numel(hB);numel(bBIIR)+numel(aBIIR)-1];
disp('COMPARACION FINAL DEL AUDIO B - NO ES SNR VERDADERA');
disp(table(Metodo,Orden_causal,Orden_efectivo,Cambio_RMS_dB, ...
    Reduccion_sobre_2kHz_dB,Correlacion_auxiliar,Mult_causal));

%% 7. GUARDADO, ESCUCHA Y CIERRE
% Nombres exactos solicitados en la asignación.
audiowrite('A_FIR.wav',normalizaSiRecorta(yA_FIR),fs);
audiowrite('A_IIR.wav',normalizaSiRecorta(yA_IIR),fs);
audiowrite('B_FIR.wav',normalizaSiRecorta(yB_FIR),fs);
audiowrite('B_IIR.wav',normalizaSiRecorta(yB_IIR),fs);

fprintf(['\nARCHIVOS GENERADOS: A_FIR.wav, A_IIR.wav, B_FIR.wav y B_IIR.wav\n' ...
 'Los FIR se implementaron causalmente y se compensó su retardo fuera de línea.\n' ...
 'Los IIR causales se analizaron, pero los WAV se generaron con filtfilt (fase cero).\n' ...
 'En B no se afirma una SNR verdadera: voz y ruido se solapan y no hay referencia limpia.\n']);

if ESCUCHAR
    fprintf('\nESCUCHA A: original -> FIR -> IIR\n'); escucharTres(xA,yA_FIR,yA_IIR,fs);
    fprintf('ESCUCHA B: original -> FIR -> IIR\n'); escucharTres(xB,yB_FIR,yB_IIR,fs);
end

%% FUNCIONES LOCALES
function x = mono(x)
if size(x,2)>1, x=mean(x,2); end
x=x(:);
end

function y = dB10(x), y=10*log10(max(abs(x),1e-20)); end
function y = dB20(x), y=20*log10(max(abs(x),1e-12)); end

function [f,Pmedio,Psegmentos] = psdSegmentada(x,fs,L)
% Hann periódica, 50 % de traslape y PSD unilateral promedio.
x=x(:)-mean(x); salto=L/2; inicios=1:salto:(numel(x)-L+1);
w=0.5-0.5*cos(2*pi*(0:L-1)'/L);
Psegmentos=zeros(L/2+1,numel(inicios));
for q=1:numel(inicios)
    s=x(inicios(q):inicios(q)+L-1).*w; X=fft(s,L);
    P=abs(X(1:L/2+1)).^2/(fs*sum(w.^2)); P(2:end-1)=2*P(2:end-1);
    Psegmentos(:,q)=P;
end
Pmedio=mean(Psegmentos,2); f=(0:L/2)'*fs/L;
end

function graficaSegmentos(f,Pseg,Pmed,limite,n,titulo)
nexttile; hold on; n=min(n,size(Pseg,2));
plot(f,dB10(Pseg(:,1:n)),'Color',[.75 .75 .75]);
plot(f,dB10(Pmed),'b','LineWidth',1.4); xlim(limite); grid on;
title(titulo); xlabel('Frecuencia [Hz]'); ylabel('PSD [dB/Hz]');
end

function [f,P] = periodogramaTramo(x,fs,L,tipo)
x=x(1:L)-mean(x(1:L));
if strcmpi(tipo,'hann'), w=0.5-0.5*cos(2*pi*(0:L-1)'/L); else, w=ones(L,1); end
X=fft(x.*w,L); P=abs(X(1:L/2+1)).^2/(fs*sum(w.^2));
P(2:end-1)=2*P(2:end-1); f=(0:L/2)'*fs/L;
end

function h = firHann(tipo,fs,frec,N)
% Método de ventanas: respuesta ideal sinc por Hann simétrica.
assert(mod(N,2)==0,'Se usa orden par para obtener retardo entero.');
m=(-N/2:N/2)'; delta=zeros(N+1,1); delta(N/2+1)=1;
switch lower(tipo)
    case 'low'
        fc=frec; hd=(2*fc/fs)*sinc((2*fc/fs)*m);
    case 'stop'
        f1=frec(1); f2=frec(2);
        hb=(2*f2/fs)*sinc((2*f2/fs)*m)-(2*f1/fs)*sinc((2*f1/fs)*m);
        hd=delta-hb;
    otherwise, error('Tipo FIR no reconocido.');
end
w=0.5-0.5*cos(2*pi*(0:N)'/N); h=hd.*w; h=(h/sum(h)).';
end

function y = firCompensado(h,x)
D=(numel(h)-1)/2; y0=filter(h,1,[x;zeros(D,1)]); y=y0(D+1:D+numel(x));
end

function c = correlacion(x,y)
N=min(numel(x),numel(y)); x=x(1:N)-mean(x(1:N)); y=y(1:N)-mean(y(1:N));
c=(x'*y)/(norm(x)*norm(y)+eps);
end

function v = cambioRMS(x,y,idx)
v=20*log10((sqrt(mean(y(idx).^2))+eps)/(sqrt(mean(x(idx).^2))+eps));
end

function A = amplitudTono(x,fs,f0)
x=x(:)-mean(x); n=(0:numel(x)-1)';
A=2/numel(x)*abs(sum(x.*exp(-1j*2*pi*f0*n/fs)));
end

function m = metricasA(x,y,fs,idx)
m=[20*log10((amplitudTono(y,fs,60)+eps)/(amplitudTono(x,fs,60)+eps)), ...
   cambioRMS(x,y,idx),correlacion(x,y)];
end

function m = metricasB(x,y,ref,fs,idx,fc)
[f,Px]=psdSegmentada(x,fs,4096); [~,Py]=psdSegmentada(y,fs,4096);
alta=f>=fc; reduccion=10*log10((sum(Py(alta))+eps)/(sum(Px(alta))+eps));
m=[cambioRMS(x,y,idx),reduccion,correlacion(ref,y)];
end

function fz = cerosCriticosFIR(h,fs,rango)
N=numel(h)-1; [H,f]=freqz(h,1,2^18,fs);
Hreal=real(H.*exp(1j*2*pi*f*(N/2)/fs));
q=f(1:end-1)>=rango(1) & f(1:end-1)<=rango(2) & Hreal(1:end-1).*Hreal(2:end)<=0;
i=find(q); fz=f(i)-Hreal(i).*(f(i+1)-f(i))./(Hreal(i+1)-Hreal(i));
end

function planoZ(z,p,limites)
th=linspace(0,2*pi,600); plot(cos(th),sin(th),'k--'); hold on;
plot(real(z),imag(z),'ob','MarkerFaceColor','none');
if ~isempty(p), plot(real(p),imag(p),'xr','LineWidth',1.2); legend('Círculo unitario','Ceros','Polos');
else, legend('Círculo unitario','Ceros'); end
axis(limites); axis equal; grid on; xlabel('Parte real'); ylabel('Parte imaginaria');
end

function b = recortaCeros(h)
tol=1e-12*max(abs(h)); q=find(abs(h)>tol); b=h(q(1):q(end));
end

function y = normaliza(x)
y=.95*x/(max(abs(x))+eps);
end

function y = normalizaSiRecorta(x)
if max(abs(x))>=1, y=normaliza(x); else, y=x; end
end

function escucharTres(x,y1,y2,fs)
sil=zeros(round(.35*fs),1); s=normaliza([x;sil;y1;sil;y2]);
sound(s,fs); pause(numel(s)/fs+.2);
end
