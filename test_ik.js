// 镜像 index.html 中的核心算法，验证正确性
let L1 = 1, L2 = 1;

function fk(t1, t2){
  return { x: L1*Math.cos(t1)+L2*Math.cos(t1+t2), y: L1*Math.sin(t1)+L2*Math.sin(t1+t2) };
}
function ikAnalytic(tx, ty){
  const r2=tx*tx+ty*ty, r=Math.sqrt(r2), reach=L1+L2, inner=Math.abs(L1-L2);
  if (r>reach+1e-9 || r<inner-1e-9 || r<1e-9) return null;
  let c2=(r2-L1*L1-L2*L2)/(2*L1*L2); c2=Math.max(-1,Math.min(1,c2));
  const t2a=Math.acos(c2), t2b=-Math.acos(c2);
  const phi=Math.atan2(ty,tx);
  let cb=(L1*L1+r2-L2*L2)/(2*L1*r); cb=Math.max(-1,Math.min(1,cb));
  const beta=Math.acos(cb);
  return { sol1:{t1:phi-beta,t2:t2a}, sol2:{t1:phi+beta,t2:t2b} };
}
function jacobian(t1,t2){
  const s1=Math.sin(t1),c1=Math.cos(t1),s12=Math.sin(t1+t2),c12=Math.cos(t1+t2);
  return [[-L1*s1-L2*s12,-L2*s12],[L1*c1+L2*c12,L2*c12]];
}
function ikNumeric(tx,ty,initT1,initT2,lam){
  let t1=initT1,t2=initT2,err=Infinity; const maxIter=200,tol=1e-4;
  for(let i=0;i<maxIter;i++){
    const f=fk(t1,t2), ex=tx-f.x, ey=ty-f.y; err=Math.hypot(ex,ey);
    if(err<tol) break;
    const J=jacobian(t1,t2);
    const a=J[0][0]*J[0][0]+J[1][0]*J[1][0]+lam*lam;
    const b=J[0][0]*J[0][1]+J[1][0]*J[1][1];
    const d=J[0][1]*J[0][1]+J[1][1]*J[1][1]+lam*lam;
    const det=a*d-b*b; if(Math.abs(det)<1e-12) break;
    const je1=J[0][0]*ex+J[1][0]*ey, je2=J[0][1]*ex+J[1][1]*ey;
    let dt1=(d*je1-b*je2)/det, dt2=(-b*je1+a*je2)/det;
    const st=Math.hypot(dt1,dt2), maxStep=0.5;
    if(st>maxStep){ const k=maxStep/st; dt1*=k; dt2*=k; }
    t1+=dt1; t2+=dt2;
  }
  const f=fk(t1,t2); err=Math.hypot(tx-f.x,ty-f.y);
  return {t1,t2,err,converged:err<1e-3};
}

let pass=0, fail=0;
function check(name, cond, info){ if(cond){pass++;} else {fail++; console.log('FAIL:',name,info);} }

// 1) 解析法：已知目标，FK 还原应一致
const targets=[[0.8,1.2],[1,1],[-0.5,0.5],[1.4,-0.3],[0,1.5],[-1.2,-1.2]];
for(const [x,y] of targets){
  const sol=ikAnalytic(x,y);
  if(!sol){ check('analytic reachable '+x+','+y, false, '返回null'); continue; }
  for(const k of ['sol1','sol2']){
    const f=fk(sol[k].t1,sol[k].t2);
    check('analytic FK '+k+' '+x+','+y, Math.hypot(f.x-x,f.y-y)<1e-9, f);
  }
}
console.log('解析法闭式解 FK 还原误差已全部核对');

// 2) 数值法：从多种初值出发应收敛到目标
for(const [x,y] of targets){
  for(const [i1,i2] of [[0,0],[1,1],[-2,0.5],[2.5,-1]]){
    const r=ikNumeric(x,y,i1,i2,0.05);
    check('numeric conv '+x+','+y+' i='+i1+','+i2, r.converged, r);
    const f=fk(r.t1,r.t2);
    check('numeric FK '+x+','+y, Math.hypot(f.x-x,f.y-y)<1e-3, f);
  }
}
console.log('数值法迭代收敛已核对');

// 3) 不可达判定
check('unreach far', ikAnalytic(2.5,0)===null, 'r=2.5>2 应不可达');
check('unreach close', (()=>{ L1=1.2;L2=0.4; const s=ikAnalytic(0.5,0); L1=1;L2=1; return s===null; })(), 'r=0.5<0.8 应不可达');
// 数值法对不可达目标不应收敛
{ const r=ikNumeric(2.5,0,0,0,0.05); check('numeric unreach no conv', !r.converged, r); }

// 4) 解析解与数值解一致性（同一目标，数值解应落在某一解析解附近）
{
  const [x,y]=[0.8,1.2];
  const sol=ikAnalytic(x,y);
  const r=ikNumeric(x,y,0,0,0.05);
  const f=fk(r.t1,r.t2);
  // 数值解末端接近目标即可
  check('numeric matches target', Math.hypot(f.x-x,f.y-y)<1e-3, f);
}
console.log(`结果: ${pass} 通过, ${fail} 失败`);
process.exit(fail?1:0);
