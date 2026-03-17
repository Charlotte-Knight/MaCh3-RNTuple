class FastTSpline3Eval {
public:
  struct Coeff {
    float x, y, b, c, d;
  };

  explicit FastTSpline3Eval(const TSpline3& spl)
  {
    const int n = spl.GetNp();
    if (n <= 0) throw std::runtime_error("TSpline3 has no points");

    coeffs_.resize(n);

    for (int i = 0; i < n; ++i) {
      double x, y, b, c, d;
      spl.GetCoeff(i, x, y, b, c, d);

      coeffs_[i] = Coeff{
        static_cast<float>(x),
        static_cast<float>(y),
        static_cast<float>(b),
        static_cast<float>(c),
        static_cast<float>(d)
      };
    }

    currSegment_ = 0;
  }

  int nPoints() const { return static_cast<int>(coeffs_.size()); }

  int findSegment(float x) const
  {
    const int n = nPoints();
    if (n <= 2) return 0;

    if (x <= coeffs_[0].x) {
      currSegment_ = 0;
      return 0;
    }

    if (x >= coeffs_[n - 1].x) {
      currSegment_ = n - 2;
      return currSegment_;
    }

    int seg = currSegment_;
    if (seg < 0) seg = 0;
    if (seg > n - 2) seg = n - 2;

    if (x >= coeffs_[seg].x && x < coeffs_[seg + 1].x) {
      currSegment_ = seg;
      return seg;
    }

    int low = 0;
    int high = n - 1;

    while (high - low > 1) {
      const int mid = (low + high) / 2;
      if (x > coeffs_[mid].x) low = mid;
      else high = mid;
    }

    seg = low;
    if (seg > n - 2) seg = n - 2;

    currSegment_ = seg;
    return seg;
  }

  float evalSegment(float x, int seg) const
  {
    const int n = nPoints();
    if (n == 0) return 0.0f;
    if (n == 1) return coeffs_[0].y;

    if (seg < 0) seg = 0;
    if (seg > n - 2) seg = n - 2;

    const Coeff& c = coeffs_[seg];
    const float dx = x - c.x;

    return fmaf(dx, fmaf(dx, fmaf(dx, c.d, c.c), c.b), c.y);
  }

  float Eval(float x) const
  {
    const int seg = findSegment(x);
    cached_value_ = evalSegment(x, seg);
    return cached_value_;
  }

  float EvalFast(float x) const
  {
    const Coeff& c = coeffs_[currSegment_];
    const float dx = x - c.x;

    cached_value_ = fmaf(dx, fmaf(dx, fmaf(dx, c.d, c.c), c.b), c.y);
    return cached_value_;
  }

  float GetCachedValue() const
  {
    return cached_value_;
  }

private:
  std::vector<Coeff> coeffs_;
  mutable int currSegment_{0};
  mutable float cached_value_{0.0f};
};