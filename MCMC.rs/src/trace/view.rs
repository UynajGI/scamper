/// Zero-copy row-major view of retained positions.
#[derive(Debug, Clone, Copy)]
pub struct TraceView<'a> {
    data: &'a [f64],
    draws: usize,
    dimension: usize,
}

impl<'a> TraceView<'a> {
    pub(crate) const fn new(data: &'a [f64], draws: usize, dimension: usize) -> Self {
        Self {
            data,
            draws,
            dimension,
        }
    }

    pub const fn draws(&self) -> usize {
        self.draws
    }

    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn draw(&self, index: usize) -> Option<&'a [f64]> {
        if index >= self.draws {
            return None;
        }
        let start = index * self.dimension;
        Some(&self.data[start..start + self.dimension])
    }

    pub fn parameter(&self, index: usize) -> Option<ParameterIter<'a>> {
        (index < self.dimension).then_some(ParameterIter {
            data: self.data,
            dimension: self.dimension,
            index,
            draw: 0,
            draws: self.draws,
        })
    }
}

/// Iterator over one parameter column of a row-major trace.
pub struct ParameterIter<'a> {
    data: &'a [f64],
    dimension: usize,
    index: usize,
    draw: usize,
    draws: usize,
}

impl Iterator for ParameterIter<'_> {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.draw >= self.draws {
            return None;
        }
        let value = self.data[self.draw * self.dimension + self.index];
        self.draw += 1;
        Some(value)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.draws - self.draw;
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for ParameterIter<'_> {}
