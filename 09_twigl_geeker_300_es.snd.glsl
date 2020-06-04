#version 430	/* version ディレクティブが必要な場合は必ず 1 行目に書くこと */
/* Copyright (C) 2020 Yosshin(@yosshin4004) */


/*
	twigl (https://twigl.app/) サウンドシェーダ互換サンプルコード。
	geeker (300 es) 対応。
*/

layout(location = 0) uniform int waveOutPosition;
#if defined(EXPORT_EXECUTABLE)
	/*
		shader minifier が compute シェーダに対応していない問題を回避するためのハック。
		以下の記述はシェーダコードとしては正しくないが、shader minifier に認識され
		minify が適用されたのち、work_around_begin: 以降のコードに置換される。
		%s は、shader minifier によるリネームが適用されたあとのシンボル名に
		置き換えらえる。
	*/
	#pragma work_around_begin:layout(std430,binding=0)buffer _{vec2 %s[];};layout(local_size_x=1)in;
	vec2 waveOutSamples[];
	#pragma work_around_end
#else
	layout(std430, binding = 0) buffer SoundOutput{ vec2 waveOutSamples[]; };
	layout(local_size_x = 1) in;
#endif


/*
	ここに twigl の mainSound 関数を張り付ける。
*/
#define bpm 132.0
#define pi2 6.2831

float calf(float i){
  return pow(2.0,i/12.0);
}
mat2 rot(float r){
  return mat2(cos(r),sin(r),-sin(r),cos(r));
}
float sigcomp1d(float w,float s){
  return (1.0/(1.0+exp(-s*w)) -0.5)*2.0;
}
float rand(vec2 st){
  return fract(sin(dot(vec2(12.9898,78.233),st))*43758.894);
}
float noise(vec2 st){
  vec2 i = floor(st);
  vec2 f = fract(st);
  float a = rand(i);
  float b = rand(i+vec2(1,0));
  float c = rand(i+vec2(0,1));
  float d = rand(i+vec2(1,1));
  vec2 u = f*f*(3.0-2.0*f);
  return mix(a,b,u.x)+(1.-u.x)*u.y*(c-a)+u.x*u.y*(d-b) -0.5;
}

float fm(float t,float f,float i,float r){
  return sin(pi2*t*f+i*sin(pi2*t*f*r));
}
float sfm(float t,float f,float i1,float r1,float i2,float r2){
  return sin(pi2*f*t+i1*fm(t,f,i2,r2)*r1);
}
float ssaw(float t,float f){
 return rand(vec2(t))*0.02+sfm(t,f,0.8,1.8,0.8,7.0)+sfm(t,f*1.008,0.8,1.8,0.8,6.0)+0.1*sfm(t,f,0.1,0.0,0.1,3.0);
}
vec2 vinegar_raw(float time,int n){
  float t = time*bpm/60.0;
  float tr = time*bpm/60.0;
  float loop = 64.0;
  float tloop = mod(t,loop);
  t = floor(2.0*t)/2.0 + pow(fract(2.0*t),1.4)/2.0;
  float hpl[32] = float[](8.,-99.,-99.,-2.,-99.,1.,-2.,-99.,
                          8.,-99.,-99.,-2.,-99.,1.,-2.,-99.,
                          -99.,-2.,-99.,1.,3.,-99.,-2.,-99.,
                          -99.,1.,-99.,-2.,3.,-99.,3.,-99.);

  float lpl[32] = float[](-99.,-99.,-2.,-99.,-99.,-99.,-99.,-2.,
                         -99.,-99.,-2.,-99.,-99.,-99.,-99.,-2.,
                        -99.,-99.,-2.,-99.,-99.,-99.,-99.,-2.,
                         -99.,-99.,-2.,-99.,-99.,-99.,-99.,-2.);

  float bl11[48] = float[](3.,3.,6.,10.,3.,6.,10.,3.,
                           3.,3.,6.,10.,3.,6.,10.,3.,
                           3.,15.,15.,10.,10.,13.,10.,13.,
                           3.,3.,15.,10.,10.,13.,10.,13.,
                           15.,15.,18.,27.,15.,18.,27.,15.,
                           15.,15.,25.,27.,15.,25.,27.,15.);

  float bl12[48] = float[](3.,-99.,6.,10.,3.,6.,10.,3.,
                           3.,-99.,6.,10.,3.,6.,10.,3.,
                           3.,-99.,15.,10.,-99.,13.,10.,13.,
                           3.,-99.,15.,10.,-99.,13.,10.,13.,
                           15.,-99.,18.,27.,15.,18.,27.,15.,
                           15.,-99.,25.,27.,15.,25.,27.,15.);

  float al[8] = float[](1.,0.,1.,1.,0.,1.,0.,1.);
  float ptoff = (n==3 || n==4)?16.:0.;
  float pt = mod(floor(t*4.0),16.0)+ptoff;
  float hpa = (hpl[int(pt)]>-50.0)?1.0:0.0;
  float lpa = (lpl[int(pt)]>-50.0)?1.0:0.0;
  float hpc=hpa*(
        fm(time,220.0*calf(hpl[int(pt)]),
           2.0*pow(fract(-4.0*t),1.5),0.79)
        +fm(time,220.0*calf(hpl[int(pt)])
           ,pow(fract(-4.0*t),15.0)*2.0,8.5))
            *fract(-t*4.0)
            *(1.0+0.3*sin(88.0*pi2*time*calf(hpl[int(pt)])));
  float lpc=lpa*(
        fm(time,110.0*calf(lpl[int(pt)]),
           2.0*pow(fract(-4.0*t),1.0),2.5)
        +fm(time,110.0*calf(lpl[int(pt)])
           ,pow(fract(-4.0*t),10.0)*3.0,7.5))
            *pow(fract(-t*4.0),0.4)
            *(1.0+0.1*sin(110.0*pi2*time*calf(lpl[int(pt)])));
  float hh1 = 0.15*(rand(vec2(time*5e2,sin(time*1e3)))
              +noise(vec2(time*5e2,sin(time*2e3))))
             *pow(sin(pi2*(t*2.0+0.61)),9.0)*(0.4+0.4*sin(pi2*t));
  float hh2 = 1.2*(noise(vec2(time*5e3,sin(time*5e3)))*0.5+rand(vec2(time)))*pow(fract(-t+0.5),25.0);

  float btoff = 0.0;
  btoff = (n==2||n==5)?16.:0.;
  btoff = (tloop>loop-4.0)?32.0:btoff;
  float bt = mod(floor(t*4.0),16.0)+btoff;
  float ba = (bl12[int(bt)]>-50.0)?1.0:0.0;
  float bacid = 1.0 - smoothstep(0.0,loop,mod(t+2.0,loop));
  float b1 = 0.0;
  float dt = 1.0/48000.0;
  float cut = 7500.0;
  float norm = 0.0;
	for(float lp=-16.0;lp<17.0;lp++){
    b1 += ssaw(time+lp*dt,55.0*calf(bl11[int(bt)]))
          *(cos(dt*lp*pi2*cut)/(1.0+pi2*dt*lp));
    norm += cos(dt*lp*pi2*cut)/(1.0+pi2*dt*lp);
  }
  b1 = (b1/norm+0.5*(1.0+bacid*pow(sin(t*pi2),2.0))*(ssaw(time,55.0*calf(bl11[int(bt)])) - b1/norm)+(1.5-0.5*bacid)*1.2*ba*(sfm(time,55.0*calf(bl12[int(bt)]),8.0*fract(-t),5.14,0.2,8.0)+noise(vec2(time*5e3,sin(time*1e4))))*pow(fract(-4.0*t),5.0)+(1.5-0.5*bacid)*ba*ssaw(time,55.0*calf(bl12[int(bt)]))*pow(fract(-4.0*t),1.0))*(1.0 + bacid*0.3*sin(t*pi2*110.0*bl11[int(bt)]));
  float bd = (0.6*(fm(tr-0.6*pow(fract(tr),2.0),62.0,3.5*pow(fract(-t),15.0),0.72))*pow(fract(-t),0.8)*log(300.0*fract(t)+1.0)+3.0*(noise(vec2(time*9e2,sin(time*1e3))))*pow(fract(-t),150.0)+0.05*(rand(vec2(time,sin(time*1e2))))*pow(fract(-t),200.0))*pow(fract(-t),1.2);
  bd =bd*0.8+0.1*sigcomp1d(bd,2.0)+0.4*sin(pi2*time*65.0)*pow(fract(-t),0.7);
  float at = mod(floor(t*4.0),8.0);
  float acid =al[int(at)]*
              fm(t+0.5*pow(sin(t),3.0),660.0,10.0*fract(-t),0.31)
              *fract(-4.0*t);
  float warp = rand(vec2(time*1e2))*0.15
               +0.8*noise(vec2(pow(fract(tr*0.25)+0.5,3.5)*2e3));

  float break1 = (tloop>loop-4.0)?0.0:1.0;
  float break2 = (tloop>loop-4.0 && tloop<loop-1.0)?1.0:0.0;
  float break3 = (tloop>loop-4.0)?1.25:1.0;

  float c1 = (n==0)?1.0:0.0;
  float c2 = (n==1)?1.0:0.0;
  float c3 = (n==2||n==5)?1.0:0.0;
  float c4 = (n==3)?1.0:0.0;
  float c5 = (n==4)?1.0:0.0;
  float c6 = (n==6)?1.0:0.0;

	//vec2 w = vec2(hpc)*0.55*rot(0.1)*break1*(c1*1.1+c2+c3+c4*1.1+c5*1.1+c6*0.6)+vec2(lpc)*0.35*rot(-0.1)*break1*(c1*1.1+c2+c3+c6*0.2)+vec2(hh1)*rot(-0.3)*break1*(c1+c2+c3+c6)+vec2(hh2)*rot(0.3)*1.2*(c3+c5+c6*0.8)+vec2(b1)*0.21*break3*(c2*1.2+c3+c4+c5)+vec2(bd)*0.6*break1*(c3+c4+c5)+vec2(acid)*0.25*rot(0.5)*break1*c3+vec2(acid*acid-0.5)*0.25*rot(-0.5)*break1*c3+vec2(warp)*break2*(c1*0.2+c2*0.8+c3+c4+c5+c6*0.4)+vec2(rand(vec2(time)))*exp(-0.3*tloop)*(c1+c6)*0.5;
	vec2 w = vec2(hpc)*0.55*rot(0.1)*break1*(c1*1.1+c2+c3+c4*1.1+c5*1.1+c6*0.6)+vec2(lpc)*0.35*rot(-0.1)*break1*(c1*1.1+c2+c3+c6*0.2);
	return w*.6;
}
vec2 mainSound(float time){
  vec2 w;
	int n = int(time*bpm/60.0/64.0);
	/*//完全に再生しない
	for(float i=-3.0;i<3.0;i++){
    w+=vinegar_raw(time+i/48000.0,n);
  }
  return vinegar_raw(time,n)*0.8 + vinegar_raw(time,n)-w/6.0;
	*/
	//途切れ途切れ再生する
	return vinegar_raw(time,n);
	//vinegar_raw内ループ回数を削減すると正常に動作
}

#define NUM_SAMPLES_PER_SEC 48000.
void main(){
	int offset = int(gl_GlobalInvocationID.x) + waveOutPosition;
	waveOutSamples[offset] = mainSound(offset / NUM_SAMPLES_PER_SEC);
}
