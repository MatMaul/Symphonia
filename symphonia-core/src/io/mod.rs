// Symphonia
// Copyright (c) 2019-2022 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! The `io` module implements composable bit- and byte-level I/O.
//!
//! The following nomenclature is used to denote where the data being read is sourced from:
//!  * A `Stream` consumes any source implementing [`ReadBytes`] one byte at a time.
//!  * A `Reader` consumes a `&[u8]`.
//!
//! The sole exception to this rule is [`MediaSourceStream`] which consumes sources implementing
//! [`MediaSource`] (aka. [`ciborium_io::Read`]).
//!
//! All `Reader`s and `Stream`s operating on bytes of data at a time implement the [`ReadBytes`]
//! trait. Likewise, all `Reader`s and `Stream`s operating on bits of data at a time implement
//! either the [`ReadBitsLtr`] or [`ReadBitsRtl`] traits depending on the order in which they
//! consume bits.

use alloc::boxed::Box;
use alloc::vec;
use core::mem;

#[cfg(feature = "std")]
pub use std::io::{Error, ErrorKind, Result, Seek, SeekFrom};

#[cfg(not(feature = "std"))]
mod no_std_io {
    use core::fmt;

    /// A minimal error kind list for no_std I/O.
    #[derive(Debug, Copy, Clone, Eq, PartialEq)]
    pub enum ErrorKind {
        /// The end of the stream was reached unexpectedly.
        UnexpectedEof,
        /// The operation was interrupted before completion.
        Interrupted,
        /// A catch-all error for other failures.
        Other,
    }

    /// A minimal no_std I/O error type.
    #[derive(Debug, Clone)]
    pub struct Error {
        kind: ErrorKind,
        message: Option<&'static str>,
    }

    impl Error {
        /// Creates a new error with a kind and static message.
        pub const fn new(kind: ErrorKind, message: &'static str) -> Self {
            Error { kind, message: Some(message) }
        }

        /// Creates a new error of kind `Other` with a static message.
        pub const fn other(message: &'static str) -> Self {
            Error::new(ErrorKind::Other, message)
        }

        /// Returns the error kind.
        pub const fn kind(&self) -> ErrorKind {
            self.kind
        }
    }

    impl From<ErrorKind> for Error {
        fn from(kind: ErrorKind) -> Self {
            Error { kind, message: None }
        }
    }

    impl fmt::Display for Error {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            if let Some(message) = self.message {
                return f.write_str(message);
            }
            let text = match self.kind {
                ErrorKind::UnexpectedEof => "unexpected end of stream",
                ErrorKind::Interrupted => "operation interrupted",
                ErrorKind::Other => "i/o error",
            };
            f.write_str(text)
        }
    }

    impl core::error::Error for Error {}

    pub type Result<T> = core::result::Result<T, Error>;
}

#[cfg(not(feature = "std"))]
pub use no_std_io::{Error, ErrorKind, Result};

use crate::io;
pub use ciborium_io::{Read, Write};

mod bit;
mod buf_reader;
mod media_source_stream;
mod monitor_stream;
mod scoped_stream;

pub use bit::*;
pub use buf_reader::BufReader;
pub use media_source_stream::{MediaSourceStream, MediaSourceStreamOptions};
pub use monitor_stream::{Monitor, MonitorStream};
pub use scoped_stream::ScopedStream;

/// `MediaSource` is a composite trait of [`ciborium_io::Read`] and, when `std` is enabled,
/// [`std::io::Seek`]. A source *must* implement this trait to be used by [`MediaSourceStream`].
///
/// Seeking is an optional capability that can be queried at runtime.
#[cfg(feature = "std")]
pub trait MediaSource: Read<Error = Error> + std::io::Read + Seek + Send + Sync {
    /// Returns if the source is seekable. This may be an expensive operation.
    fn is_seekable(&self) -> bool;

    /// Returns the length in bytes, if available. This may be an expensive operation.
    fn byte_len(&self) -> Option<u64>;
}

#[cfg(not(feature = "std"))]
pub trait MediaSource: Read<Error = Error> + Send + Sync {
    /// Returns if the source is seekable. This may be an expensive operation.
    fn is_seekable(&self) -> bool;

    /// Returns the length in bytes, if available. This may be an expensive operation.
    fn byte_len(&self) -> Option<u64>;
}

#[cfg(feature = "std")]
impl MediaSource for std::fs::File {
    /// Returns if the `std::io::File` backing the `MediaSource` is seekable.
    ///
    /// Note: This operation involves querying the underlying file descriptor for information and
    /// may be moderately expensive. Therefore it is recommended to cache this value if used often.
    fn is_seekable(&self) -> bool {
        // If the file's metadata is available, and the file is a regular file (i.e., not a FIFO,
        // etc.), then the MediaSource will be seekable. Otherwise assume it is not. Note that
        // metadata() follows symlinks.
        match self.metadata() {
            Ok(metadata) => metadata.is_file(),
            _ => false,
        }
    }

    /// Returns the length in bytes of the `std::io::File` backing the `MediaSource`.
    ///
    /// Note: This operation involves querying the underlying file descriptor for information and
    /// may be moderately expensive. Therefore it is recommended to cache this value if used often.
    fn byte_len(&self) -> Option<u64> {
        match self.metadata() {
            Ok(metadata) => Some(metadata.len()),
            _ => None,
        }
    }
}

#[cfg(feature = "std")]
impl<T: std::convert::AsRef<[u8]> + Send + Sync> MediaSource for std::io::Cursor<T> {
    /// Always returns true since a `std::io::Cursor<u8>` is always seekable.
    fn is_seekable(&self) -> bool {
        true
    }

    /// Returns the length in bytes of the `std::io::Cursor<u8>` backing the `MediaSource`.
    fn byte_len(&self) -> Option<u64> {
        // Get the underlying container, usually &Vec<T>.
        let inner = self.get_ref();
        // Get slice from the underlying container, &[T], for the len() function.
        Some(inner.as_ref().len() as u64)
    }
}

/// `ReadOnlySource` wraps any source implementing [`ciborium_io::Read`] in an unseekable
/// [`MediaSource`].
pub struct ReadOnlySource<R: Read<Error = Error>> {
    inner: R,
}

impl<R: Read<Error = Error> + Send> ReadOnlySource<R> {
    /// Instantiates a new `ReadOnlySource<R>` by taking ownership and wrapping the provided
    /// `Read`er.
    pub fn new(inner: R) -> Self {
        ReadOnlySource { inner }
    }

    /// Gets a reference to the underlying reader.
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Gets a mutable reference to the underlying reader.
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Unwraps this `ReadOnlySource<R>`, returning the underlying reader.
    pub fn into_inner(self) -> R {
        self.inner
    }
}

#[cfg(feature = "std")]
impl<R: Read<Error = Error> + std::io::Read + Send + Sync> MediaSource for ReadOnlySource<R> {
    fn is_seekable(&self) -> bool {
        false
    }

    fn byte_len(&self) -> Option<u64> {
        None
    }
}

#[cfg(not(feature = "std"))]
impl<R: Read<Error = Error> + Send + Sync> MediaSource for ReadOnlySource<R> {
    fn is_seekable(&self) -> bool {
        false
    }

    fn byte_len(&self) -> Option<u64> {
        None
    }
}

#[cfg(not(feature = "std"))]
impl<R: Read<Error = Error>> Read for ReadOnlySource<R> {
    type Error = Error;

    fn read_exact(&mut self, buf: &mut [u8]) -> Result<()> {
        self.inner.read_exact(buf)
    }
}

#[cfg(feature = "std")]
impl<R: Read<Error = Error> + std::io::Read> std::io::Read for ReadOnlySource<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        std::io::Read::read(&mut self.inner, buf)
    }
}

#[cfg(feature = "std")]
impl<R: Read<Error = Error> + std::io::Read> Seek for ReadOnlySource<R> {
    fn seek(&mut self, _: SeekFrom) -> Result<u64> {
        Err(Error::other("source does not support seeking"))
    }
}

/// `ReadBytes` provides methods to read bytes and interpret them as little- or big-endian
/// unsigned integers or floating-point values of standard widths.
pub trait ReadBytes {
    /// Reads a single byte from the stream and returns it or an error.
    fn read_byte(&mut self) -> io::Result<u8>;

    /// Reads two bytes from the stream and returns them in read-order or an error.
    fn read_double_bytes(&mut self) -> io::Result<[u8; 2]>;

    /// Reads three bytes from the stream and returns them in read-order or an error.
    fn read_triple_bytes(&mut self) -> io::Result<[u8; 3]>;

    /// Reads four bytes from the stream and returns them in read-order or an error.
    fn read_quad_bytes(&mut self) -> io::Result<[u8; 4]>;

    /// Reads up-to the number of bytes required to fill buf or returns an error.
    fn read_buf(&mut self, buf: &mut [u8]) -> io::Result<usize>;

    /// Reads exactly the number of bytes required to fill be provided buffer or returns an error.
    fn read_buf_exact(&mut self, buf: &mut [u8]) -> io::Result<()>;

    /// Reads a single unsigned byte from the stream and returns it or an error.
    #[inline(always)]
    fn read_u8(&mut self) -> io::Result<u8> {
        self.read_byte()
    }

    /// Reads a single signed byte from the stream and returns it or an error.
    #[inline(always)]
    fn read_i8(&mut self) -> io::Result<i8> {
        Ok(self.read_byte()? as i8)
    }

    /// Reads two bytes from the stream and interprets them as an unsigned 16-bit little-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_u16(&mut self) -> io::Result<u16> {
        Ok(u16::from_le_bytes(self.read_double_bytes()?))
    }

    /// Reads two bytes from the stream and interprets them as an signed 16-bit little-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_i16(&mut self) -> io::Result<i16> {
        Ok(i16::from_le_bytes(self.read_double_bytes()?))
    }

    /// Reads two bytes from the stream and interprets them as an unsigned 16-bit big-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_be_u16(&mut self) -> io::Result<u16> {
        Ok(u16::from_be_bytes(self.read_double_bytes()?))
    }

    /// Reads two bytes from the stream and interprets them as an signed 16-bit big-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_be_i16(&mut self) -> io::Result<i16> {
        Ok(i16::from_be_bytes(self.read_double_bytes()?))
    }

    /// Reads three bytes from the stream and interprets them as an unsigned 24-bit little-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_u24(&mut self) -> io::Result<u32> {
        let mut buf = [0u8; mem::size_of::<u32>()];
        buf[0..3].clone_from_slice(&self.read_triple_bytes()?);
        Ok(u32::from_le_bytes(buf))
    }

    /// Reads three bytes from the stream and interprets them as an signed 24-bit little-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_i24(&mut self) -> io::Result<i32> {
        Ok(((self.read_u24()? << 8) as i32) >> 8)
    }

    /// Reads three bytes from the stream and interprets them as an unsigned 24-bit big-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_be_u24(&mut self) -> io::Result<u32> {
        let mut buf = [0u8; mem::size_of::<u32>()];
        buf[0..3].clone_from_slice(&self.read_triple_bytes()?);
        Ok(u32::from_be_bytes(buf) >> 8)
    }

    /// Reads three bytes from the stream and interprets them as an signed 24-bit big-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_be_i24(&mut self) -> io::Result<i32> {
        Ok(((self.read_be_u24()? << 8) as i32) >> 8)
    }

    /// Reads four bytes from the stream and interprets them as an unsigned 32-bit little-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_u32(&mut self) -> io::Result<u32> {
        Ok(u32::from_le_bytes(self.read_quad_bytes()?))
    }

    /// Reads four bytes from the stream and interprets them as an signed 32-bit little-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_i32(&mut self) -> io::Result<i32> {
        Ok(i32::from_le_bytes(self.read_quad_bytes()?))
    }

    /// Reads four bytes from the stream and interprets them as an unsigned 32-bit big-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_be_u32(&mut self) -> io::Result<u32> {
        Ok(u32::from_be_bytes(self.read_quad_bytes()?))
    }

    /// Reads four bytes from the stream and interprets them as a signed 32-bit big-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_be_i32(&mut self) -> io::Result<i32> {
        Ok(i32::from_be_bytes(self.read_quad_bytes()?))
    }

    /// Reads eight bytes from the stream and interprets them as an unsigned 64-bit little-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_u64(&mut self) -> io::Result<u64> {
        let mut buf = [0u8; mem::size_of::<u64>()];
        self.read_buf_exact(&mut buf)?;
        Ok(u64::from_le_bytes(buf))
    }

    /// Reads eight bytes from the stream and interprets them as an signed 64-bit little-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_i64(&mut self) -> io::Result<i64> {
        let mut buf = [0u8; mem::size_of::<i64>()];
        self.read_buf_exact(&mut buf)?;
        Ok(i64::from_le_bytes(buf))
    }

    /// Reads eight bytes from the stream and interprets them as an unsigned 64-bit big-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_be_u64(&mut self) -> io::Result<u64> {
        let mut buf = [0u8; mem::size_of::<u64>()];
        self.read_buf_exact(&mut buf)?;
        Ok(u64::from_be_bytes(buf))
    }

    /// Reads eight bytes from the stream and interprets them as an signed 64-bit big-endian
    /// integer or returns an error.
    #[inline(always)]
    fn read_be_i64(&mut self) -> io::Result<i64> {
        let mut buf = [0u8; mem::size_of::<i64>()];
        self.read_buf_exact(&mut buf)?;
        Ok(i64::from_be_bytes(buf))
    }

    /// Reads four bytes from the stream and interprets them as a 32-bit little-endian IEEE-754
    /// floating-point value.
    #[inline(always)]
    fn read_f32(&mut self) -> io::Result<f32> {
        Ok(f32::from_le_bytes(self.read_quad_bytes()?))
    }

    /// Reads four bytes from the stream and interprets them as a 32-bit big-endian IEEE-754
    /// floating-point value.
    #[inline(always)]
    fn read_be_f32(&mut self) -> io::Result<f32> {
        Ok(f32::from_be_bytes(self.read_quad_bytes()?))
    }

    /// Reads four bytes from the stream and interprets them as a 64-bit little-endian IEEE-754
    /// floating-point value.
    #[inline(always)]
    fn read_f64(&mut self) -> io::Result<f64> {
        let mut buf = [0u8; mem::size_of::<u64>()];
        self.read_buf_exact(&mut buf)?;
        Ok(f64::from_le_bytes(buf))
    }

    /// Reads four bytes from the stream and interprets them as a 64-bit big-endian IEEE-754
    /// floating-point value.
    #[inline(always)]
    fn read_be_f64(&mut self) -> io::Result<f64> {
        let mut buf = [0u8; mem::size_of::<u64>()];
        self.read_buf_exact(&mut buf)?;
        Ok(f64::from_be_bytes(buf))
    }

    /// Reads up-to the number of bytes requested, and returns a boxed slice of the data or an
    /// error.
    fn read_boxed_slice(&mut self, len: usize) -> io::Result<Box<[u8]>> {
        let mut buf = vec![0u8; len];
        let actual_len = self.read_buf(&mut buf)?;
        buf.truncate(actual_len);
        Ok(buf.into_boxed_slice())
    }

    /// Reads exactly the number of bytes requested, and returns a boxed slice of the data or an
    /// error.
    fn read_boxed_slice_exact(&mut self, len: usize) -> io::Result<Box<[u8]>> {
        let mut buf = vec![0u8; len];
        self.read_buf_exact(&mut buf)?;
        Ok(buf.into_boxed_slice())
    }

    /// Reads bytes from the stream into a supplied buffer until a byte pattern is matched. Returns
    /// a mutable slice to the valid region of the provided buffer.
    #[inline(always)]
    fn scan_bytes<'a>(&mut self, pattern: &[u8], buf: &'a mut [u8]) -> io::Result<&'a mut [u8]> {
        self.scan_bytes_aligned(pattern, 1, buf)
    }

    /// Reads bytes from a stream into a supplied buffer until a byte patter is matched on an
    /// aligned byte boundary. Returns a mutable slice to the valid region of the provided buffer.
    fn scan_bytes_aligned<'a>(
        &mut self,
        pattern: &[u8],
        align: usize,
        buf: &'a mut [u8],
    ) -> io::Result<&'a mut [u8]>;

    /// Ignores the specified number of bytes from the stream or returns an error.
    fn ignore_bytes(&mut self, count: u64) -> io::Result<()>;

    /// Gets the position of the stream.
    fn pos(&self) -> u64;
}

impl<R: ReadBytes> ReadBytes for &mut R {
    #[inline(always)]
    fn read_byte(&mut self) -> io::Result<u8> {
        (*self).read_byte()
    }

    #[inline(always)]
    fn read_double_bytes(&mut self) -> io::Result<[u8; 2]> {
        (*self).read_double_bytes()
    }

    #[inline(always)]
    fn read_triple_bytes(&mut self) -> io::Result<[u8; 3]> {
        (*self).read_triple_bytes()
    }

    #[inline(always)]
    fn read_quad_bytes(&mut self) -> io::Result<[u8; 4]> {
        (*self).read_quad_bytes()
    }

    #[inline(always)]
    fn read_buf(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        (*self).read_buf(buf)
    }

    #[inline(always)]
    fn read_buf_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        (*self).read_buf_exact(buf)
    }

    #[inline(always)]
    fn scan_bytes_aligned<'a>(
        &mut self,
        pattern: &[u8],
        align: usize,
        buf: &'a mut [u8],
    ) -> io::Result<&'a mut [u8]> {
        (*self).scan_bytes_aligned(pattern, align, buf)
    }

    #[inline(always)]
    fn ignore_bytes(&mut self, count: u64) -> io::Result<()> {
        (*self).ignore_bytes(count)
    }

    #[inline(always)]
    fn pos(&self) -> u64 {
        (**self).pos()
    }
}

impl<S: SeekBuffered> SeekBuffered for &mut S {
    fn ensure_seekback_buffer(&mut self, len: usize) {
        (*self).ensure_seekback_buffer(len)
    }

    fn unread_buffer_len(&self) -> usize {
        (**self).unread_buffer_len()
    }

    fn read_buffer_len(&self) -> usize {
        (**self).read_buffer_len()
    }

    fn seek_buffered(&mut self, pos: u64) -> u64 {
        (*self).seek_buffered(pos)
    }

    fn seek_buffered_rel(&mut self, delta: isize) -> u64 {
        (*self).seek_buffered_rel(delta)
    }
}

/// `SeekBuffered` provides methods to seek within the buffered portion of a stream.
pub trait SeekBuffered {
    /// Ensures that `len` bytes will be available for backwards seeking if `len` bytes have been
    /// previously read.
    fn ensure_seekback_buffer(&mut self, len: usize);

    /// Get the number of bytes buffered but not yet read.
    ///
    /// Note: This is the maximum number of bytes that can be seeked forwards within the buffer.
    fn unread_buffer_len(&self) -> usize;

    /// Gets the number of bytes buffered and read.
    ///
    /// Note: This is the maximum number of bytes that can be seeked backwards within the buffer.
    fn read_buffer_len(&self) -> usize;

    /// Seek within the buffered data to an absolute position in the stream. Returns the position
    /// seeked to.
    fn seek_buffered(&mut self, pos: u64) -> u64;

    /// Seek within the buffered data relative to the current position in the stream. Returns the
    /// position seeked to.
    ///
    /// The range of `delta` is clamped to the inclusive range defined by
    /// `-read_buffer_len()..=unread_buffer_len()`.
    fn seek_buffered_rel(&mut self, delta: isize) -> u64;

    /// Seek backwards within the buffered data.
    ///
    /// This function is identical to [`SeekBuffered::seek_buffered_rel`] when a negative delta is
    /// provided.
    fn seek_buffered_rev(&mut self, delta: usize) {
        assert!(delta < isize::MAX as usize);
        self.seek_buffered_rel(-(delta as isize));
    }
}

impl<F: FiniteStream> FiniteStream for &mut F {
    fn byte_len(&self) -> u64 {
        (**self).byte_len()
    }

    fn bytes_read(&self) -> u64 {
        (**self).bytes_read()
    }

    fn bytes_available(&self) -> u64 {
        (**self).bytes_available()
    }
}

/// A `FiniteStream` is a stream that has a known length in bytes.
pub trait FiniteStream {
    /// Returns the length of the the stream in bytes.
    fn byte_len(&self) -> u64;

    /// Returns the number of bytes that have been read.
    fn bytes_read(&self) -> u64;

    /// Returns the number of bytes available for reading.
    fn bytes_available(&self) -> u64;
}
