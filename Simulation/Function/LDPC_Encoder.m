function encodedData = LDPC_Encoder(data, g)
    data = gf(data,2);
    g = gf(g,2);
    encodedData = data*g;
    encodedData = (encodedData.x);
end