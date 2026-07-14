document.addEventListener('DOMContentLoaded', function () {
  const sections = document.querySelectorAll('.variable-section');
  sections.forEach(function (section) {
    const tabs = section.querySelectorAll('.tab');
    const panels = section.querySelectorAll('.panel');
    tabs.forEach(function (tab) {
      tab.addEventListener('click', function () {
        const target = tab.getAttribute('data-projection');
        // Update this section's tabs only — sections are independent,
        // so switching projection on temperature doesn't touch wind, etc.
        tabs.forEach(function (t) {
          const isMatch = t === tab;
          t.classList.toggle('is-active', isMatch);
          t.setAttribute('aria-selected', isMatch ? 'true' : 'false');
        });
        panels.forEach(function (panel) {
          const isMatch = panel.getAttribute('data-projection') === target;
          panel.classList.toggle('is-active', isMatch);
          // Pause any playing video in panels that just became hidden,
          // so switching tabs doesn't leave it running in the background.
          if (!isMatch) {
            const video = panel.querySelector('video');
            if (video && !video.paused) video.pause();
          }
        });
      });
    });
  });

  // Playback speed controls — each panel with a video gets its own
  // Slow/Normal/Fast buttons, scoped so clicking one only affects
  // the video in that same panel.
  const speedButtons = document.querySelectorAll('.speed-btn');
  speedButtons.forEach(function (btn) {
    btn.addEventListener('click', function () {
      const panel = btn.closest('.panel');
      const video = panel ? panel.querySelector('video') : null;
      if (!video) return;

      video.playbackRate = parseFloat(btn.getAttribute('data-speed'));

      // Only toggle active state among this panel's own speed buttons.
      const siblingButtons = panel.querySelectorAll('.speed-btn');
      siblingButtons.forEach(function (b) {
        b.classList.toggle('is-active', b === btn);
      });
    });
  });
});
