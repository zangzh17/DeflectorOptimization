function LayerTextures = TexturesGen1D(PatternIn,NLevel,dx,Nz,nTop,nDevice)
% Define textures for each layer
LayerTextures = cell(1,Nz+2);
LayerTextures{1} = {nTop};
LayerTextures{2} = {nBot};
dz = (NLevel-1)/Nz;
for ii=1:Nz
    h = (ii-1) * dz;
    index = (PatternIn(:)'>h);
    index = [index,~index(end)];
    index = diff(index);
    lowIndexPos = find(index==1)*dx;
    highIndexPos = find(index==-1)*dx;
    index = [ones(size(lowIndexPos))*nTop, ones(size(highIndexPos))*nDevice];
    pos = [lowIndexPos,highIndexPos];
    [pos, I] = sort(pos); index = index(I);
    LayerTextures{ii+2} = {pos, index};
end
end