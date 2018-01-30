function X = circshift_ad_2D(X, k)

    X = circshift_ad(circshift_ad(X, k(1))', k(2))';

end
